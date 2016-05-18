// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2016 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//
/*
 *
 *  Created on: Jun 3, 2015
 *      Author: jduddy
 */

#include "indel_overlapper.hh"
#include "blt_util/log.hh"
#include "ScoringModelManager.hh"

//#define DEBUG_GVCF



indel_overlapper::indel_overlapper(const ScoringModelManager& model, const reference_contig_segment& ref, std::shared_ptr<variant_pipe_stage_base> destination)
    : variant_pipe_stage_base(destination)
    , _CM(model)
    , _ref(ref)
    , _indel_end_pos(0)
{
    // this component doesn't make any sense without a destination:
    assert(destination);
}



void indel_overlapper::process(std::unique_ptr<site_info> site)
{
    std::unique_ptr<digt_site_info> si(downcast<digt_site_info>(std::move(site)));

    // resolve any current or previous indels before queuing site:
    if (si->pos>=_indel_end_pos)
    {
        process_overlaps();
    }
    else
    {
        _site_buffer.push_back(std::move(si));
        return;
    }

    assert(si->pos>=_indel_end_pos);
    assert(_nonvariant_indel_buffer.empty());

    _sink->process(std::move(si));
}



static
bool
is_het_indel(const starling_diploid_indel_core& dindel)
{
    return (dindel.max_gt==STAR_DIINDEL::HET);
}

static
bool
check_is_nonvariant_indel(const starling_diploid_indel_core& dindel)
{
    return (dindel.max_gt==STAR_DIINDEL::NOINDEL);
}

void indel_overlapper::process(std::unique_ptr<indel_info> indel)
{
    auto ii(downcast<digt_indel_info>(std::move(indel)));

    auto& call(ii->first());

    // we can't handle breakends at all right now:
    if (call._ik.is_breakpoint()) return;

    const bool is_nonvariant_indel = check_is_nonvariant_indel(call._dindel);

    // don't handle homozygous reference calls unless genotyping is forced
    if (is_nonvariant_indel && !call._dindel.is_forced_output) return;

    if (ii->pos>_indel_end_pos)
    {
        process_overlaps();
    }

    if (is_nonvariant_indel)
    {
        _nonvariant_indel_buffer.push_back(std::move(ii));
    }
    else
    {
        _indel_end_pos=std::max(_indel_end_pos,call._ik.right_pos());
        _indel_buffer.push_back(std::move(ii));
    }
}

static
bool
is_simple_indel_overlap(
    const reference_contig_segment& ref,
    const std::vector<std::unique_ptr<digt_indel_info>>& indel_buffer)
{
#ifdef DEBUG_GVCF
    log_os << "CHIRP: " << __FUNCTION__ << " START\n";
#endif

    // check for very very simple overlap condition -- these are the cases that are easy to
    // glue together, although many more non-simple cases could be resolved if we wanted to
    // put in the work
    //
    // check for 2 overlapping hets:

    if (indel_buffer.size() != 2) return false;

    const digt_indel_info& ii0(*indel_buffer[0]);
    const digt_indel_info& ii1(*indel_buffer[1]);

    const digt_indel_call& ic0(ii0.first());
    const digt_indel_call& ic1(ii1.first());

    const bool isTwoHets = (is_het_indel(ic0._dindel) && is_het_indel(ic1._dindel));

#ifdef DEBUG_GVCF
    log_os << "CHIRP: " << __FUNCTION__ << " isTwoHets: " << isTwoHets << "\n";
#endif

    if (! isTwoHets) return false;

    // also make sure we don't create two matching alt sequences
    // (two matching ALTs will fail vcf output

    // there's going to be 1 (possibly empty) fill range in front of one haplotype
    // and one possibly empty fill range on the back of one haplotype
    const pos_t indel_end_pos=std::max(ic1._ik.right_pos(),ic0._ik.right_pos());
    const pos_t indel_begin_pos(ii0.pos-1);

    // get the VCF ALT string associated with overlapping indel:
    auto get_overlap_alt = [&] (const digt_indel_info& ii)
    {
        std::string leading_seq,trailing_seq;
        const auto& ic(ii.first());
        // extend leading sequence start back 1 for vcf compat, and end back 1 to concat with vcf_indel_seq
        ref.get_substring(indel_begin_pos,(ii.pos-indel_begin_pos)-1,leading_seq);
        const unsigned trail_len(indel_end_pos-ic._ik.right_pos());
        ref.get_substring(indel_end_pos-trail_len,trail_len,trailing_seq);

        return leading_seq + ic._iri.vcf_indel_seq + trailing_seq;
    };
    const std::string alt0 = get_overlap_alt(*indel_buffer[0]);
    const std::string alt1 = get_overlap_alt(*indel_buffer[1]);

#ifdef DEBUG_GVCF
    log_os << "CHIRP: " << __FUNCTION__ << " alt0: " << alt0 << " alt1: " << alt1 << "\n";
#endif

    return (alt0 != alt1);
}



void
indel_overlapper::
dump(std::ostream& os) const
{
    os << "indel_overlapper:"
       << " nSites: " << _site_buffer.size()
       << " nIndels: " << _indel_buffer.size()
       << " indel_end_pos: " << _indel_end_pos << "\n";
    os << "buffered sites:\n";
    for (const auto& site : _site_buffer)
    {
        os << *site << "\n";
    }

    os << "buffered indels:\n";
    for (const auto& indel : _indel_buffer)
    {
        indel->dump(os);
        os << "\n";
    }
}



void
indel_overlapper::
process_overlaps()
{
    try
    {
        process_overlaps_impl();
    }
    catch (...)
    {
        log_os << "ERROR: exception caught in process_overlaps()\n";
        dump(log_os);
        throw;
    }
}


namespace VARQUEUE
{
enum index_t
{
    NONE,
    INDEL,
    NONVARIANT_INDEL,
    SITE
};
}


// this doesn't really generalize or tidy up the (implicit) indel/site priority queue, but
// just dumps the ugliness into one place:
static
VARQUEUE::index_t
nextVariantType(
    const std::vector<std::unique_ptr<digt_indel_info>>& indel_buffer,
    const std::vector<std::unique_ptr<digt_indel_info>>& nonvariant_indel_buffer,
    const std::vector<std::unique_ptr<digt_site_info>>& site_buffer,
    const unsigned indel_index,
    const unsigned nonvariant_indel_index,
    const unsigned site_index)
{
    const bool is_indel(indel_index<indel_buffer.size());
    const bool is_nonvariant_indel(nonvariant_indel_index<nonvariant_indel_buffer.size());
    const bool is_site(site_index<site_buffer.size());

    if ((!is_indel) && (!is_nonvariant_indel) && (!is_site))
    {
        return VARQUEUE::NONE;
    }

    const bool AlessB(is_indel && ((! is_nonvariant_indel) || (indel_buffer[indel_index]->pos <= nonvariant_indel_buffer[nonvariant_indel_index]->pos)));
    const bool AlessC(is_indel && ((! is_site) || (indel_buffer[indel_index]->pos <= site_buffer[site_index]->pos)));
    const bool BlessC(is_nonvariant_indel && ((! is_site) || (nonvariant_indel_buffer[nonvariant_indel_index]->pos <= site_buffer[site_index]->pos)));

    if (AlessB)
    {
        if (AlessC)
        {
            return VARQUEUE::INDEL;
        }
        else
        {
            return VARQUEUE::SITE;
        }
    }
    else
    {
        if (BlessC)
        {
            return VARQUEUE::NONVARIANT_INDEL;
        }
        else
        {
            return VARQUEUE::SITE;
        }
    }
}



void indel_overlapper::process_overlaps_impl()
{
#ifdef DEBUG_GVCF
    log_os << "CHIRP: " << __FUNCTION__ << " START\n";
#endif

    if (_indel_buffer.empty() && _nonvariant_indel_buffer.empty()) return;

    bool is_conflict(false);

    // do standard or overlap indel processing:
    if (_indel_buffer.size()==1)
    {
        // simple case of no overlap:
        modify_single_indel_record(*_indel_buffer[0]);
    }
    else if (_indel_buffer.size() > 1)
    {
        if (is_simple_indel_overlap(_ref,_indel_buffer))
        {
            // handle the simplest possible overlap case (two hets):
            modify_overlap_indel_record();
        }
        else
        {
            // mark the whole region as conflicting
            modify_conflict_indel_record();
            is_conflict=true;
        }
    }

    // simple processing for all nonvariant_indels:
    for (auto& nonvariant_indel : _nonvariant_indel_buffer)
    {
        modify_single_indel_record(*nonvariant_indel);
    }

    // process sites to be consistent with overlapping indels:
    for (auto& si : _site_buffer)
    {
#ifdef DEBUG_GVCF
        log_os << "CHIRP: indel overlapping site: " << si->pos << "\n";
#endif
        modify_overlapping_site(*_indel_buffer[0], *si, _CM);
    }

    unsigned indel_index(0);
    unsigned nonvariant_indel_index(0);
    unsigned site_index(0);

    // order all buffered indel and site record output according to VCF formatting rules:
    while (true)
    {
        const VARQUEUE::index_t nextvar = nextVariantType(_indel_buffer,_nonvariant_indel_buffer,_site_buffer,indel_index,nonvariant_indel_index,site_index);

        if      (nextvar == VARQUEUE::NONE)
        {
            break;
        }
        else if (nextvar == VARQUEUE::INDEL)
        {
            _sink->process(std::move(_indel_buffer[indel_index]));
            if (is_conflict)
            {
                // emit each conflict record
                indel_index++;
            }
            else
            {
                // just emit the overlapped or single non-conflict record
                indel_index=_indel_buffer.size();
            }
        }
        else if (nextvar == VARQUEUE::NONVARIANT_INDEL)
        {
            _sink->process(std::move(_nonvariant_indel_buffer[nonvariant_indel_index]));
            nonvariant_indel_index++;
        }
        else if (nextvar == VARQUEUE::SITE)
        {
            _sink->process(std::move(_site_buffer[site_index]));
            site_index++;
        }
        else
        {
            assert(false && "unexpected varqueue type");
        }
    }

    _indel_buffer.clear();
    _nonvariant_indel_buffer.clear();
    _site_buffer.clear();
}



void indel_overlapper::modify_overlapping_site(const digt_indel_info& ii, digt_site_info& si, const ScoringModelManager& model)
{
    const pos_t offset(si.pos-ii.pos);
    assert(offset>=0);

    if (ii.first().filters.test(VCF_FILTERS::IndelConflict))
    {
        modify_indel_conflict_site(si);
    }
    else
    {
        modify_indel_overlap_site( ii,
                                   ii.get_ploidy(offset),
                                   si, model);
    }
}

// set the CIGAR string:
void
indel_overlapper::
modify_single_indel_record(digt_indel_info& ii)
{
    ii.first().set_hap_cigar();
    _CM.classify_indel(ii, ii.first());
}



void indel_overlapper::modify_indel_overlap_site(
    const digt_indel_info& ii,
    const unsigned ploidy,
    digt_site_info& si,
    const ScoringModelManager& model)
{
#ifdef DEBUG_GVCF
    log_os << "CHIRP: indel_overlap_site smod before: " << si.smod << "\n";
#endif

    // if overlapping indel has any filters, mark as site conflict
    // (note that we formerly had the site inherit indel filters, but
    // this interacts poorly with empirical scoring)
    if (! ii.first().filters.none())
    {
        si.smod.set_filter(VCF_FILTERS::SiteConflict);
    }

#ifdef DEBUG_GVCF
    log_os << "CHIRP: indel_overlap_site smod after: " << si.smod << "\n";
#endif

    // limit qual and gq values to those of the indel
    si.dgt.genome.snp_qphred = std::min(si.dgt.genome.snp_qphred,ii.first()._dindel.indel_qphred);
    si.smod.gqx = std::min(si.smod.gqx,ii.first()._dindel.max_gt_qphred);

    // change ploidy:
    if (ploidy==1)
    {
        if (DIGT::is_het(si.smod.max_gt))
        {
            si.smod.set_filter(VCF_FILTERS::SiteConflict);
            //si.smod.modified_gt=MODIFIED_SITE_GT::UNKNOWN;
        }
        else
        {
            if (si.smod.max_gt == si.dgt.ref_gt)
            {
                si.smod.modified_gt=MODIFIED_SITE_GT::ZERO;
            }
            else
            {
                si.smod.modified_gt=MODIFIED_SITE_GT::ONE;
            }
        }
    }
    else if (ploidy==0)
    {
        if (si.smod.max_gt == si.dgt.ref_gt)
        {
            si.smod.modified_gt=MODIFIED_SITE_GT::UNKNOWN;
            si.smod.is_zero_ploidy=true;
            if (si.dgt.is_noploid())
            {
                si.smod.unset_filter(VCF_FILTERS::PloidyConflict);
            }
        }
        else
        {
            si.smod.set_filter(VCF_FILTERS::SiteConflict);
        }
    }
    else if (ploidy!=2)
    {
        assert(false && "Unexpected ploidy value");
    }

    // after all those changes we need to rerun the site filters:
    model.classify_site(si, si.smod);
}



void indel_overlapper::modify_indel_conflict_site(digt_site_info& si)
{
    si.smod.set_filter(VCF_FILTERS::IndelConflict);
}

// figure out the per-site ploidy inside of indel based on each haplotype's match descriptor:





void
indel_overlapper::modify_overlap_indel_record()
{
#ifdef DEBUG_GVCF
    log_os << "CHIRP: " << __FUNCTION__ << " START\n";
#endif

    // can only handle simple 2-indel overlaps right now:
    assert(_indel_buffer.size()==2);

    // TODO: hackey
    for (auto& ii : _indel_buffer)
        ii->_is_overlap = true;

    _CM.classify_indels(_indel_buffer);

    _indel_buffer[0]->add_overlap(_ref, *_indel_buffer[1]);
}



// set the CIGAR string:
void
indel_overlapper::modify_conflict_indel_record()
{
#ifdef DEBUG_GVCF
    log_os << "CHIRP: " << __FUNCTION__ << " START\n";
#endif

    assert(_indel_buffer.size()>1);

    for (auto& ii : _indel_buffer)
    {
        ii->first().set_hap_cigar();

        ii->first().set_filter(VCF_FILTERS::IndelConflict);

        _CM.classify_indel(*ii, ii->first());
    }
}



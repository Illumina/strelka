// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Starka
// Copyright (c) 2009-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//
/*
 *
 *  Created on: Jun 3, 2015
 *      Author: jduddy
 */

#include "indel_overlapper.hh"
#include "calibration_models.hh"

//#define DEBUG_GVCF


#ifdef DEBUG_GVCF
#include "blt_util/log.hh"
#endif





indel_overlapper::indel_overlapper(const calibration_models& model, const reference_contig_segment& ref, std::shared_ptr<variant_pipe_stage_base> destination)
    : variant_pipe_stage_base(destination)
    , _CM(model)
    , _ref(ref)
    , _indel_end_pos(0)
{

}

void indel_overlapper::flush()
{
    // flush out accumulated sites & indels
    process_overlaps();
    variant_pipe_stage_base::flush();
}

void indel_overlapper::process(std::unique_ptr<site_info> site)
{
    std::unique_ptr<digt_site_info> si(downcast<digt_site_info>(std::move(site)));

    // resolve any current or previous indels before queuing site:
    if (! _indel_buffer.empty())
    {
        if (si->pos>=_indel_end_pos)
        {
            process_overlaps();
        }
        else
        {
            _site_buffer.push_back(std::move(si));
            return;
        }
    }
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
is_no_indel(const starling_diploid_indel_core& dindel)
{
    return (dindel.max_gt==STAR_DIINDEL::NOINDEL);
}

void indel_overlapper::process(std::unique_ptr<indel_info> indel)
{
    auto ii(downcast<digt_indel_info>(std::move(indel)));

    auto& call(ii->first());

    // we can't handle breakends at all right now:
    if (call._ik.is_breakpoint()) return;

    // don't handle homozygous reference calls unless genotyping is forced
    if (is_no_indel(call._dindel) && !call._dindel.is_forced_output) return;


    bool no_indel = is_no_indel(call._dindel);

    if ((! _indel_buffer.empty()) && ((ii->pos>_indel_end_pos) || no_indel))
    {
        process_overlaps();
    }
    _indel_end_pos=std::max(_indel_end_pos,call._ik.right_pos());
    _indel_buffer.push_back(std::move(ii));
    // clear the current homRef indel
    if (no_indel)
    {
        process_overlaps();
    }
}

static
bool
is_simple_indel_overlap(
    const reference_contig_segment& ref,
    const std::vector<std::unique_ptr<digt_indel_info>>& indel_buffer)
{
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

    return (alt0 != alt1);
}




void indel_overlapper::process_overlaps()
{
    if (0==_indel_buffer.size()) return;

    bool is_conflict(false);

    // do the overlap processing:
    if (_indel_buffer.size()==1)
    {
        // simple case of no overlap:
        modify_single_indel_record();
    }
    else
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

    //    *_osptr << "INDEL_SIZE: " << _indel_buffer.size() << "\n";

    // process sites to be consistent with overlapping indels:
    for (auto& si : _site_buffer)
    {
#ifdef DEBUG_GVCF
        log_os << "CHIRP: indel overlapping site: " << si->pos << "\n";
#endif
        modify_overlapping_site(*_indel_buffer[0], *si, _CM);
    }

    unsigned indel_index(0);
    unsigned site_index(0);

    while (true)
    {
        const bool is_indel(indel_index<_indel_buffer.size());
        const bool is_site(site_index<_site_buffer.size());
        if (! (is_indel || is_site)) break;

        if (is_indel && ((! is_site) || _indel_buffer[indel_index]->pos <= _site_buffer[site_index]->pos))
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
        else
        {
            // emit site:
            //log_os << "site record" << "\n";
            _sink->process(std::move(_site_buffer[site_index]));
            site_index++;
        }
    }
    _indel_buffer.clear();
    _site_buffer.clear();
}

void indel_overlapper::modify_overlapping_site(const digt_indel_info& ii, digt_site_info& si, const calibration_models& model)
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
indel_overlapper::modify_single_indel_record()
{
    assert(_indel_buffer.size()==1);

    digt_indel_info& ii(*_indel_buffer[0]);
    ii.first().set_hap_cigar();

    _CM.clasify_indel(ii, ii.first());
}



void indel_overlapper::modify_indel_overlap_site(
    const digt_indel_info& ii,
    const unsigned ploidy,
    digt_site_info& si,
    const calibration_models& model)
{
#ifdef DEBUG_GVCF
    log_os << "CHIRP: indel_overlap_site smod before: " << si.smod << "\n";
    log_os << "CHIRP: indel_overlap_site imod before: " << ii.imod << "\n";
#endif

    // if overlapping indel has any filters, mark as site conflict
    // (note that we formerly had the site inherit indel filters, but
    // this interacts poorly with VQSR)
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
    model.clasify_site(si, si.smod);
}



void indel_overlapper::modify_indel_conflict_site(digt_site_info& si)
{
    si.smod.set_filter(VCF_FILTERS::IndelConflict);
}

// figure out the per-site ploidy inside of indel based on each haplotype's match descriptor:





void
indel_overlapper::modify_overlap_indel_record()
{
    // can only handle simple 2-indel overlaps right now:
    assert(_indel_buffer.size()==2);

    // TODO: hackey
    for (auto& ii : _indel_buffer)
        ii->_is_overlap = true;

    _CM.clasify_indels(_indel_buffer);

    _indel_buffer[0]->add_overlap(_ref, *_indel_buffer[1]);

}



// set the CIGAR string:
void
indel_overlapper::modify_conflict_indel_record()
{
    assert(_indel_buffer.size()>1);

    for (auto& ii : _indel_buffer)
    {
        ii->first().set_hap_cigar();

        ii->first().set_filter(VCF_FILTERS::IndelConflict);

        _CM.clasify_indel(*ii, ii->first());
    }
}



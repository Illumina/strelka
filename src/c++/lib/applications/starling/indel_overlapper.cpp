// -*- mode: c++; indent-tabs-mode: nil; -*-
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





indel_overlapper::indel_overlapper(const calibration_models& model, const reference_contig_segment& ref, pos_t head_pos, variant_pipe_stage& destination)
    : variant_pipe_stage(destination)
    , _CM(model)
    , _ref(ref)
    , _indel_end_pos(0)
    , _head_pos(head_pos)
{

}

void indel_overlapper::flush()
{
    // flush out accumulated sites & indels
    process_overlaps();
    variant_pipe_stage::flush();
}

void indel_overlapper::process(site_info& si)
{
    if (si.smod.is_phased_region)
    {
        _head_pos=si.pos+si.phased_ref.length();
    }
    else
    {
        _head_pos=si.pos+1;
    }

    // resolve any current or previous indels before queuing site:
    if (! _indel_buffer.empty())
    {
        if (si.pos>=_indel_end_pos)
        {
            process_overlaps();
        }
        else
        {
            _site_buffer.push_back(si);
            return;
        }
    }
    _sink->process(si);

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

void indel_overlapper::process(indel_info& ii)
{
    // we can't handle breakends at all right now:
    if (ii.ik().is_breakpoint()) return;

    // don't handle homozygous reference calls unless genotyping is forced
    if (is_no_indel(ii.dindel) && !ii.dindel.is_forced_output) return;


    if ((! _indel_buffer.empty()) && ((ii.pos>_indel_end_pos) || is_no_indel(ii.dindel)))
    {
        process_overlaps();
    }
    _indel_buffer.push_back(ii);
    _indel_end_pos=std::max(_indel_end_pos,ii.ik().right_pos());
     // clear the current homRef indel
     if (is_no_indel(ii.dindel))
     {
         process_overlaps();
     }
}

static
bool
is_simple_indel_overlap(const std::vector<indel_info>& indel_buffer)
{
    return (indel_buffer.size()==2 &&
            is_het_indel(indel_buffer[0].dindel) &&
            is_het_indel(indel_buffer[1].dindel));
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
        if (is_simple_indel_overlap(_indel_buffer))
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
    for (site_info& si : _site_buffer)
    {
#ifdef DEBUG_GVCF
        log_os << "CHIRP: indel overlapping site: " << it->pos << "\n";
#endif
        const pos_t offset(si.pos-_indel_buffer[0].pos);
        assert(offset>=0);
        if (! is_conflict)
        {
            modify_indel_overlap_site( _indel_buffer[0],
                                       _indel_buffer[0].get_ploidy(offset),
                                       si);
        }
        else
        {
            modify_indel_conflict_site(si);
        }
    }

    unsigned indel_index(0);
    unsigned site_index(0);

    // TODO: redo the data structures such that overlapped indels can be passed as a single record to process()
    while (true)
    {
        const bool is_indel(indel_index<_indel_buffer.size());
        const bool is_site(site_index<_site_buffer.size());
        if (! (is_indel || is_site)) break;

        if (is_indel && ((! is_site) || _indel_buffer[indel_index].pos <= _site_buffer[site_index].pos))
        {
            _sink->process(_indel_buffer[indel_index]);
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
            _sink->process(_site_buffer[site_index]);
            site_index++;
        }
    }
    _indel_buffer.clear();
    _site_buffer.clear();
}

// set the CIGAR string:
void
indel_overlapper::modify_single_indel_record()
{
    assert(_indel_buffer.size()==1);

    indel_info& ii(_indel_buffer[0]);
    ii.set_hap_cigar();

    _CM.clasify_indel(ii, ii.imod());
}



void indel_overlapper::modify_indel_overlap_site(
    const indel_info& ii,
    const unsigned ploidy,
    site_info& si)
{
#ifdef DEBUG_GVCF
    log_os << "CHIRP: indel_overlap_site smod before: " << si.smod << "\n";
    log_os << "CHIRP: indel_overlap_site imod before: " << ii.imod << "\n";
#endif

    // if overlapping indel has any filters, mark as site conflict
    // (note that we formerly had the site inherit indel filters, but
    // this interacts poorly with VQSR)
    if (! ii.imod().filters.none())
    {
        si.smod.set_filter(VCF_FILTERS::SiteConflict);
    }

#ifdef DEBUG_GVCF
    log_os << "CHIRP: indel_overlap_site smod after: " << si.smod << "\n";
#endif

    // limit qual and gq values to those of the indel
    si.dgt.genome.snp_qphred = std::min(si.dgt.genome.snp_qphred,ii.dindel.indel_qphred);
    si.smod.gqx = std::min(si.smod.gqx,ii.dindel.max_gt_qphred);

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
    _CM.clasify_site(si, si.smod);
}



void indel_overlapper::modify_indel_conflict_site(site_info& si)
{
    si.smod.set_filter(VCF_FILTERS::IndelConflict);
}

// figure out the per-site ploidy inside of indel based on each haplotype's match descriptor:





void
indel_overlapper::modify_overlap_indel_record()
{
    // can only handle simple 2-indel overlaps right now:
    assert(_indel_buffer.size()==2);

    _indel_buffer[0].add_overlap(_ref, _indel_buffer[1]);

    // TODO: classification uses the overlap stuff - fix setting that
    // by moving the merging to the output phase
    _CM.clasify_indels(_indel_buffer);
}



// set the CIGAR string:
void
indel_overlapper::modify_conflict_indel_record()
{
    assert(_indel_buffer.size()>1);

    for (auto it = _indel_buffer.begin(); it != _indel_buffer.end(); ++it)
    {
        indel_info& ii(*it);
        ii.set_hap_cigar();

        ii._imod.front().set_filter(VCF_FILTERS::IndelConflict);

        _CM.clasify_indel(ii, ii.imod());
    }
}



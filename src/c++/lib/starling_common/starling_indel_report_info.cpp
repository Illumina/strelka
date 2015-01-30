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

/// \file
///
/// \author Chris Saunders
///


#include "blt_common/ref_context.hh"
#include "blt_util/seq_util.hh"
#include "starling_common/pos_basecall_buffer.hh"
#include "starling_common/starling_indel_report_info.hh"

#include "boost/lexical_cast.hpp"

#include <cassert>
#include <cmath>


#ifdef DEBUG_REPORT
#include "blt_util/log.hh"
#endif



///
static
void
get_indel_desc(const indel_key& ik,
               std::string& indel_desc)
{

    indel_desc.clear();

    if ((ik.type == INDEL::INSERT) ||
        (ik.type == INDEL::DELETE) ||
        (ik.type == INDEL::SWAP))
    {
        indel_desc += boost::lexical_cast<std::string>(ik.length);
        if ((ik.type == INDEL::INSERT) ||
            (ik.type == INDEL::SWAP))
        {
            indel_desc += 'I';
            if (ik.swap_dlength>0)
            {
                indel_desc += boost::lexical_cast<std::string>(ik.swap_dlength);
                indel_desc += 'D';
            }
        }
        else
        {
            indel_desc += 'D';
        }
    }
    else if (ik.type == INDEL::BP_LEFT)
    {
        indel_desc = "BP_LEFT";
    }
    else if (ik.type == INDEL::BP_RIGHT)
    {
        indel_desc = "BP_RIGHT";
    }
    else
    {
        assert(0);
    }
}


#if 0
static
void
copy_subseq(const std::string& in_seq,
            const pos_t start_pos,
            const pos_t end_pos,
            std::string& out_seq)
{

    const char* ip(in_seq.c_str());
    const pos_t is(in_seq.size());
    out_seq.clear();
    for (pos_t p(start_pos); p<end_pos; ++p)
    {
        out_seq += get_seq_base(ip,is,p);
    }
}
#endif



static
void
copy_ref_subseq(const reference_contig_segment& ref,
                const pos_t start_pos,
                const pos_t end_pos,
                std::string& out_seq)
{

    out_seq.clear();
    for (pos_t p(start_pos); p<end_pos; ++p)
    {
        out_seq += ref.get_base(p);
    }
}



static
void
set_delete_seq(const indel_key& ik,
               const reference_contig_segment& ref,
               std::string& seq)
{
    copy_ref_subseq(ref,ik.pos,ik.right_pos(),seq);
}



/// get the indel cigar and ref and indel strings used in the indel
/// summary line output
///
static
void
get_indel_summary_strings(const indel_key& ik,
                          const indel_data& id,
                          const reference_contig_segment& ref,
                          std::string& indel_desc,
                          std::string& indel_seq,
                          std::string& indel_ref_seq)
{

    get_indel_desc(ik,indel_desc);

    if ((ik.type == INDEL::INSERT) ||
        (ik.type == INDEL::BP_LEFT) ||
        (ik.type == INDEL::BP_RIGHT))
    {
        indel_seq = id.get_insert_seq();
        indel_ref_seq = std::string(indel_seq.size(),'-');
    }
    else if (ik.type == INDEL::DELETE)
    {
        set_delete_seq(ik,ref,indel_ref_seq);
        indel_seq = std::string(indel_ref_seq.size(),'-');
    }
    else if (ik.type == INDEL::SWAP)
    {
        indel_seq = id.get_insert_seq();
        set_delete_seq(ik,ref,indel_ref_seq);
        int idelt(static_cast<int>(ik.length)-static_cast<int>(ik.swap_dlength));
        if (idelt>0)
        {
            indel_ref_seq += std::string(idelt,'-');
        }
        else
        {
            indel_seq += std::string(-idelt,'-');
        }
    }
    else
    {
        assert(0);
    }

    assert(indel_seq.size() == indel_ref_seq.size());
}



/// get the indel cigar and ref and indel strings used in the indel
/// summary line output
///
static
void
get_vcf_summary_strings(const indel_key& ik,
                        const indel_data& id,
                        const reference_contig_segment& ref,
                        std::string& vcf_indel_seq,
                        std::string& vcf_ref_seq)
{

    if       (ik.is_breakpoint())
    {
        if       (ik.type == INDEL::BP_LEFT)
        {
            copy_ref_subseq(ref,ik.pos-1,ik.pos,vcf_ref_seq);
            vcf_indel_seq = vcf_ref_seq + id.get_insert_seq() + '.';
        }
        else if (ik.type == INDEL::BP_RIGHT)
        {
            copy_ref_subseq(ref,ik.pos,ik.pos+1,vcf_ref_seq);
            vcf_indel_seq = '.' + id.get_insert_seq() + vcf_ref_seq;
        }
        else
        {
            assert(0);
        }
    }
    else
    {
        copy_ref_subseq(ref,ik.pos-1,ik.pos+ik.delete_length(),vcf_ref_seq);
        copy_ref_subseq(ref,ik.pos-1,ik.pos,vcf_indel_seq);
        vcf_indel_seq += id.get_insert_seq();
    }
}



static
void
set_repeat_info(const indel_key& ik,
                const reference_contig_segment& ref,
                starling_indel_report_info& iri)
{
    if (! ((iri.it == INDEL::INSERT) ||
           (iri.it == INDEL::DELETE) ||
           (iri.it == INDEL::SWAP))) return;

    unsigned insert_repeat_count(0);
    unsigned delete_repeat_count(0);

    if       (iri.it == INDEL::INSERT)
    {
        get_seq_repeat_unit(iri.indel_seq,iri.repeat_unit,insert_repeat_count);
    }
    else if (iri.it == INDEL::DELETE)
    {
        get_seq_repeat_unit(iri.ref_seq,iri.repeat_unit,delete_repeat_count);
    }
    else if (iri.it == INDEL::SWAP)
    {
        std::string insert_ru;
        std::string delete_ru;
        get_seq_repeat_unit(iri.indel_seq,insert_ru,insert_repeat_count);
        get_seq_repeat_unit(iri.ref_seq,delete_ru,delete_repeat_count);
        if ((insert_ru != delete_ru) || insert_ru.empty()) return;

        iri.repeat_unit=insert_ru;
    }
    else
    {
        assert(false && "Unexpected indel type");
    }

    // count repeats in contextual sequence:
    unsigned indel_context_repeat_count(0);
    {
        const pos_t indel_begin_pos(ik.pos);
        const pos_t indel_end_pos(ik.right_pos());
        const int repeat_unit_size(static_cast<int>(iri.repeat_unit.size()));

        // count upstream repeats:
        for (pos_t i(indel_begin_pos-repeat_unit_size); i>=0; i-=repeat_unit_size)
        {
            bool is_repeat(true);
            for (int j(0); j<repeat_unit_size; ++j)
            {
                if (ref.get_base(i+j) != iri.repeat_unit[j])
                {
                    is_repeat = false;
                    break;
                }
            }
            if (! is_repeat) break;
            indel_context_repeat_count += 1;
        }

        // count downstream repeats:
        const pos_t rs(ref.end());
        for (pos_t i(indel_end_pos); (i+static_cast<pos_t>(repeat_unit_size)-1)<rs; i+=repeat_unit_size)
        {
            bool is_repeat(true);
            for (int j(0); j<repeat_unit_size; ++j)
            {
                if (ref.get_base(i+j) != iri.repeat_unit[j])
                {
                    is_repeat = false;
                    break;
                }
            }
            if (! is_repeat) break;
            indel_context_repeat_count += 1;
        }
    }

    iri.ref_repeat_count = indel_context_repeat_count+delete_repeat_count;
    iri.indel_repeat_count = indel_context_repeat_count+insert_repeat_count;
}



void
get_starling_indel_report_info(const indel_key& ik,
                               const indel_data& id,
                               const reference_contig_segment& ref,
                               starling_indel_report_info& iri)
{
    // indel summary info
    get_indel_summary_strings(ik,id,ref,iri.desc,iri.indel_seq,iri.ref_seq);
    get_vcf_summary_strings(ik,id,ref,iri.vcf_indel_seq,iri.vcf_ref_seq);

    iri.it=ik.type;

    const pos_t indel_begin_pos(ik.pos);
    const pos_t indel_end_pos(ik.right_pos());

    // reference context:
    {
        static const unsigned INDEL_CONTEXT_SIZE(10);

        if (ik.type != INDEL::BP_RIGHT)
        {
            iri.ref_upstream.clear();
            for (pos_t i(indel_begin_pos-static_cast<pos_t>(INDEL_CONTEXT_SIZE)); i<indel_begin_pos; ++i)
            {
                iri.ref_upstream += ref.get_base(i);
            }
        }
        else
        {
            iri.ref_upstream = "N/A";
        }
        if (ik.type != INDEL::BP_LEFT)
        {
            iri.ref_downstream.clear();
            for (pos_t i(indel_end_pos); i<(indel_end_pos+static_cast<pos_t>(INDEL_CONTEXT_SIZE)); ++i)
            {
                iri.ref_downstream += ref.get_base(i);
            }
        }
        else
        {
            iri.ref_downstream = "N/A";
        }
    }

    // repeat analysis:
    set_repeat_info(ik,ref,iri);

    // interrupted hpol compuation:
    iri.ihpol=get_interrupted_hpol_size(indel_begin_pos-1,ref);
    iri.ihpol=std::max(iri.ihpol,get_interrupted_hpol_size(indel_begin_pos,ref));
    if (indel_begin_pos != indel_end_pos)
    {
        iri.ihpol=std::max(iri.ihpol,get_interrupted_hpol_size(indel_end_pos-1,ref));
        iri.ihpol=std::max(iri.ihpol,get_interrupted_hpol_size(indel_end_pos,ref));
    }
}



read_path_scores
indel_lnp_to_pprob(const starling_base_deriv_options& dopt,
                   const read_path_scores& path_lnp,
                   const bool is_tier2_pass,
                   const bool is_use_alt_indel)
{
    unsigned n_alleles(2);
    if (is_use_alt_indel)
    {
        n_alleles += path_lnp.alt_indel.size();
    }

    static const double allele_prior(1./static_cast<double>(n_alleles));
    static const double allele_lnprior(std::log(allele_prior));

    read_path_scores pprob;
    read_path_scores::score_t pprob_nonsite = dopt.get_nonsite_path_lnp(is_tier2_pass,path_lnp.nsite) + dopt.nonsite_lnprior;
    pprob.ref     = path_lnp.ref     + dopt.site_lnprior + allele_lnprior;
    pprob.indel   = path_lnp.indel   + dopt.site_lnprior + allele_lnprior;

    if (is_use_alt_indel)
    {
        for (const auto& val : path_lnp.alt_indel)
        {
            pprob.alt_indel.push_back(std::make_pair(val.first,(val.second + dopt.site_lnprior + allele_lnprior)));
        }
    }

    double scale(std::max(pprob_nonsite,std::max(pprob.ref,pprob.indel)));

    if (is_use_alt_indel)
    {
        for (const auto& val : pprob.alt_indel)
        {
            if (scale < val.second) scale = val.second;
        }
    }

    pprob_nonsite = std::exp(pprob_nonsite-scale);
    pprob.ref = std::exp(pprob.ref-scale);
    pprob.indel = std::exp(pprob.indel-scale);
    if (is_use_alt_indel)
    {
        for (auto& val : pprob.alt_indel)
        {
            val.second = std::exp((val.second)-scale);
        }
    }


    double sum(pprob_nonsite+pprob.ref+pprob.indel);
    if (is_use_alt_indel)
    {
        for (const auto& val : pprob.alt_indel)
        {
            sum += val.second;
        }
    }

    // pprob_nonsite /= sum; /// no point in normalizing this if we aren't adding it back into pprob
    pprob.ref /= sum;
    pprob.indel /= sum;
    if (is_use_alt_indel)
    {
        for (auto& val : pprob.alt_indel)
        {
            val.second /= sum;
        }
    }

    return pprob;
}



void
get_starling_indel_sample_report_info(const starling_base_deriv_options& dopt,
                                      const indel_key& ik,
                                      const indel_data& id,
                                      const pos_basecall_buffer& bc_buff,
                                      const bool is_tier2_pass,
                                      const bool is_use_alt_indel,
                                      starling_indel_sample_report_info& isri)
{
    // get read info:
    {
        static const double path_pprob_thresh(0.999);

        unsigned n_subscore_reads(0);

        typedef indel_data::score_t::const_iterator siter;
        siter i(id.read_path_lnp.begin()), i_end(id.read_path_lnp.end());
        for (; i!=i_end; ++i)
        {
            const read_path_scores& path_lnp(i->second);

            // optionally skip tier2 data:
            if ((! is_tier2_pass) && (! path_lnp.is_tier1_read)) continue;

            const read_path_scores pprob(indel_lnp_to_pprob(dopt,path_lnp,is_tier2_pass,is_use_alt_indel));
            if       (pprob.ref >= path_pprob_thresh)
            {
                isri.n_q30_ref_reads++;
            }
            else if (pprob.indel >= path_pprob_thresh)
            {
                isri.n_q30_indel_reads++;
            }
            else
            {
                typedef read_path_scores::alt_indel_t::const_iterator aciter;

                bool is_alt_found(false);
#if 0
                if (pprob.is_alt && (pprob.alt >= path_pprob_thresh))
                {
                    isri.n_q30_alt_reads++;
                    is_alt_found=true;
                }
#else
                aciter j(pprob.alt_indel.begin()), j_end(pprob.alt_indel.end());
                for (; j!=j_end; ++j)
                {
                    if (j->second >= path_pprob_thresh)
                    {
                        isri.n_q30_alt_reads++;
                        is_alt_found=true;
                        break;
                    }
                }
#endif
                if (! is_alt_found)
                {
                    n_subscore_reads++;
                }
            }
        }

        // total number of reads with non-zero, yet insufficient indel
        // breakpoint overlap
        const unsigned n_suboverlap_tier1_reads(id.suboverlap_tier1_read_ids.size());
        isri.n_other_reads = (n_subscore_reads+n_suboverlap_tier1_reads);

        if (is_tier2_pass)
        {
            const unsigned n_suboverlap_tier2_reads(id.suboverlap_tier2_read_ids.size());
            isri.n_other_reads += n_suboverlap_tier2_reads;
        }
    }

    {
        // get depth of indel:
        pos_t depth_pos(ik.pos-1);
        if (ik.type==INDEL::BP_RIGHT) depth_pos=ik.pos;
        const snp_pos_info& spi(bc_buff.get_pos(depth_pos));
        isri.depth=spi.calls.size();
    }
}

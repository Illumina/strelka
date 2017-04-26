//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2017 Illumina, Inc.
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

///
/// \author Chris Saunders
///


#include "AlleleReportInfoUtil.hh"
#include "blt_common/ref_context.hh"
#include "blt_util/seq_util.hh"
#include "starling_common/pos_basecall_buffer.hh"
#if 0
#include "starling_common/AlleleReportInfo.hh"

#include "boost/lexical_cast.hpp"

#include <cassert>
#include <cmath>
#endif


#ifdef DEBUG_REPORT
#include "blt_util/log.hh"
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



void
getSingleIndelAlleleVcfSummaryStrings(
    const IndelKey& indelKey,
    const IndelData& indelData,
    const reference_contig_segment& ref,
    std::string& vcf_indel_seq,
    std::string& vcf_ref_seq)
{
    if       (indelKey.is_breakpoint())
    {
        if       (indelKey.type == INDEL::BP_LEFT)
        {
            copy_ref_subseq(ref,indelKey.pos-1,indelKey.pos,vcf_ref_seq);
            vcf_indel_seq = vcf_ref_seq + indelData.getBreakpointInsertSeq() + '.';
        }
        else if (indelKey.type == INDEL::BP_RIGHT)
        {
            copy_ref_subseq(ref,indelKey.pos,indelKey.pos+1,vcf_ref_seq);
            vcf_indel_seq = '.' + indelData.getBreakpointInsertSeq() + vcf_ref_seq;
        }
        else
        {
            assert(0);
        }
    }
    else
    {
        copy_ref_subseq(ref,indelKey.pos-1,indelKey.pos+indelKey.delete_length(),vcf_ref_seq);
        copy_ref_subseq(ref,indelKey.pos-1,indelKey.pos,vcf_indel_seq);
        vcf_indel_seq += indelKey.insert_seq();
    }
}

// Find shortest prefix of unencoded sequence that has not been encountered in the encoded sequence: (TODO: this can be made more efficient)
unsigned
shortest_unencountered(
    const reference_contig_segment& ref,
    const pos_t pos,
    const unsigned num_encoded,
    const bool left = false)
{
    std::string encoded;
    std::string new_prefix;
    // get substring encoded thus far:
    if (left)
    {
        ref.get_substring(pos-num_encoded+1,num_encoded,encoded);
    }
    else
    {
        ref.get_substring(pos,num_encoded,encoded);
    }
    unsigned len(0);
    do
    {
        ++len;
        if (left)
        {
            ref.get_substring(pos-num_encoded-len+1,len,new_prefix);
        }
        else
        {
            ref.get_substring(pos+num_encoded,len,new_prefix);
        }
    }
    // keep looking while string of size len has already been encountered:
    while (encoded.find(new_prefix) != std::string::npos);
    return len;
}

unsigned
compute_context_compressability(
    const reference_contig_segment& ref,
    const pos_t left_pos,
    const pos_t right_pos,
    const unsigned numkeys)
{
    // complexity to the left:
    unsigned num_encoded_left(1);
    for (unsigned key(0); key<numkeys-1; ++key)
    {
        num_encoded_left += shortest_unencountered(ref,left_pos-1,num_encoded_left,true);
    }

    // complexity to the right:
    unsigned num_encoded_right(1);
    for (unsigned key(0); key<numkeys-1; ++key)
    {
        num_encoded_right += shortest_unencountered(ref,right_pos,num_encoded_right,false);
    }

    return std::max(num_encoded_left, num_encoded_right);
}

static
void
set_repeat_info(
    const IndelKey& indelKey,
    const reference_contig_segment& ref,
    AlleleReportInfo& indelReportInfo)
{
    // Encoded lengths using fixed number of Zev-Lempel 1977 keywords. See http://www.lptmc.jussieu.fr/user/lesne/PRE-Short.pdf :
    indelReportInfo.contextCompressability = compute_context_compressability(ref, indelKey.pos, indelKey.right_pos(), 5);

    if (! ((indelReportInfo.it == SimplifiedIndelReportType::INSERT) ||
           (indelReportInfo.it == SimplifiedIndelReportType::DELETE) ||
           (indelReportInfo.it == SimplifiedIndelReportType::SWAP))) return;

    unsigned insert_repeat_count(0);
    unsigned delete_repeat_count(0);

    if      (indelReportInfo.it == SimplifiedIndelReportType::INSERT)
    {
        get_seq_repeat_unit(indelKey.insert_seq(), indelReportInfo.repeat_unit, insert_repeat_count);
    }
    else if (indelReportInfo.it == SimplifiedIndelReportType::DELETE)
    {
        std::string deletedSeq;
        copy_ref_subseq(ref,indelKey.pos, indelKey.pos+indelKey.delete_length(), deletedSeq);

        get_seq_repeat_unit(deletedSeq, indelReportInfo.repeat_unit, delete_repeat_count);
    }
    else if (indelReportInfo.it == SimplifiedIndelReportType::SWAP)
    {
        std::string insert_ru;
        get_seq_repeat_unit(indelKey.insert_seq(), insert_ru, insert_repeat_count);

        std::string delete_ru;
        std::string deletedSeq;
        copy_ref_subseq(ref,indelKey.pos, indelKey.pos+indelKey.delete_length(), deletedSeq);
        get_seq_repeat_unit(deletedSeq, delete_ru, delete_repeat_count);

        if ((insert_ru != delete_ru) || insert_ru.empty()) return;

        indelReportInfo.repeat_unit=insert_ru;
    }
    else
    {
        assert(false && "Unexpected indel type");
    }
    indelReportInfo.repeat_unit_length=indelReportInfo.repeat_unit.size();

    // count repeats in contextual sequence:
    unsigned indel_context_repeat_count(0);
    {
        const pos_t indel_begin_pos(indelKey.pos);
        const pos_t indel_end_pos(indelKey.right_pos());
        const int repeat_unit_size(static_cast<int>(indelReportInfo.repeat_unit.size()));

        // count upstream repeats:
        for (pos_t i(indel_begin_pos-repeat_unit_size); i>=0; i-=repeat_unit_size)
        {
            bool is_repeat(true);
            for (int j(0); j<repeat_unit_size; ++j)
            {
                if (ref.get_base(i+j) != indelReportInfo.repeat_unit[j])
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
                if (ref.get_base(i+j) != indelReportInfo.repeat_unit[j])
                {
                    is_repeat = false;
                    break;
                }
            }
            if (! is_repeat) break;
            indel_context_repeat_count += 1;
        }
    }

    indelReportInfo.ref_repeat_count = indel_context_repeat_count+delete_repeat_count;
    indelReportInfo.indel_repeat_count = indel_context_repeat_count+insert_repeat_count;
}



void
getAlleleReportInfo(
    const IndelKey& indelKey,
    const reference_contig_segment& ref,
    AlleleReportInfo& indelReportInfo)
{
    indelReportInfo.it=SimplifiedIndelReportType::getRateType(indelKey);

    const pos_t indel_begin_pos(indelKey.pos);
    const pos_t indel_end_pos(indelKey.right_pos());

    // repeat analysis:
    set_repeat_info(indelKey,ref,indelReportInfo);

    // interrupted hpol compuation:
    indelReportInfo.ihpol=get_interrupted_hpol_size(indel_begin_pos-1,ref);
    indelReportInfo.ihpol=std::max(indelReportInfo.ihpol,get_interrupted_hpol_size(indel_begin_pos,ref));
    if (indel_begin_pos != indel_end_pos)
    {
        indelReportInfo.ihpol=std::max(indelReportInfo.ihpol,get_interrupted_hpol_size(indel_end_pos-1,ref));
        indelReportInfo.ihpol=std::max(indelReportInfo.ihpol,get_interrupted_hpol_size(indel_end_pos,ref));
    }
}



ReadPathScores
indel_lnp_to_pprob(
    const starling_base_deriv_options& dopt,
    const ReadPathScores& path_lnp,
    const bool is_tier2_pass,
    const bool is_use_alt_indel)
{
    unsigned n_alleles(2);
    if (is_use_alt_indel)
    {
        n_alleles += path_lnp.alt_indel.size();
    }

    const double allele_prior(1./static_cast<double>(n_alleles));
    const double allele_lnprior(std::log(allele_prior));

    ReadPathScores pprob;
    ReadPathScores::score_t pprob_nonsite = dopt.get_nonsite_path_lnp(is_tier2_pass,path_lnp.nsite) + dopt.nonsite_lnprior;
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
getAlleleSampleReportInfo(
    const starling_base_options& opt,
    const starling_base_deriv_options& dopt,
    const IndelKey& indelKey,
    const IndelSampleData& indelSampleData,
    const pos_basecall_buffer& bc_buff,
    const bool is_tier2_pass,
    const bool is_use_alt_indel,
    AlleleSampleReportInfo& isri)
{
    // get read info:
    {
        unsigned n_subscore_reads(0);

        for (const auto& val : indelSampleData.read_path_lnp)
        {
            const ReadPathScores& path_lnp(val.second);

            // optionally skip tier2 data:
            if ((! is_tier2_pass) && (! path_lnp.is_tier1_read)) continue;

            const ReadPathScores pprob(indel_lnp_to_pprob(dopt,path_lnp,is_tier2_pass,is_use_alt_indel));

            const unsigned usableReadPos(std::max(0,static_cast<int>(path_lnp.read_pos)));
            if       (pprob.ref >= opt.readConfidentSupportThreshold.numval())
            {
                isri.n_confident_ref_reads++;
                if (path_lnp.is_fwd_strand)
                {
                    ++isri.n_confident_ref_reads_fwd;
                }
                else
                {
                    ++isri.n_confident_ref_reads_rev;
                }

                isri.readpos_ranksum.add_observation(true, usableReadPos);
            }
            else if (pprob.indel >= opt.readConfidentSupportThreshold.numval())
            {
                isri.n_confident_indel_reads++;
                if (path_lnp.is_fwd_strand)
                {
                    ++isri.n_confident_indel_reads_fwd;
                }
                else
                {
                    ++isri.n_confident_indel_reads_rev;
                }

                isri.readpos_ranksum.add_observation(false, usableReadPos);

                {
                    // this may help the mean edge distance better generalize over sequencing assay read length:
                    static const int16_t maxDistanceFromEdge(20);

                    // add this read edge feature only for reads supporting the indel
                    isri.distanceFromReadEdge.addObservation(std::min(maxDistanceFromEdge, path_lnp.distanceFromClosestReadEdge));
                }
            }
            else
            {
                bool is_alt_found(false);
#if 0
                if (pprob.is_alt && (pprob.alt >= path_confident_support_threshold))
                {
                    isri.n_confident_alt_reads++;
                    is_alt_found=true;
                }
#else
                for (const auto& palt : pprob.alt_indel)
                {
                    if (palt.second >= opt.readConfidentSupportThreshold.numval())
                    {
                        isri.n_confident_alt_reads++;
                        if (path_lnp.is_fwd_strand)
                        {
                            ++isri.n_confident_alt_reads_fwd;
                        }
                        else
                        {
                            ++isri.n_confident_alt_reads_rev;
                        }
                        is_alt_found=true;
                        break;
                    }
                }
#endif
                if (! is_alt_found)
                {
                    n_subscore_reads++;
                    if (path_lnp.is_fwd_strand)
                    {
                        ++isri.n_other_reads_fwd;
                    }
                    else
                    {
                        ++isri.n_other_reads_rev;
                    }
                }
            }
        }

        // total number of reads with non-zero, yet insufficient indel
        // breakpoint overlap
        const unsigned n_suboverlap_tier1_reads(indelSampleData.suboverlap_tier1_read_ids.size());
        isri.n_other_reads = (n_subscore_reads+n_suboverlap_tier1_reads);

        if (is_tier2_pass)
        {
            const unsigned n_suboverlap_tier2_reads(indelSampleData.suboverlap_tier2_read_ids.size());
            isri.n_other_reads += n_suboverlap_tier2_reads;
        }
    }

    {
        // get various indel stats from the pileup:
        pos_t depth_pos(indelKey.pos-1);
        if (indelKey.type==INDEL::BP_RIGHT) depth_pos=indelKey.pos;
        const snp_pos_info& spi(bc_buff.get_pos(depth_pos));
        isri.tier1Depth=spi.calls.size();
        isri.mapqTracker=spi.mapqTracker;
    }
}

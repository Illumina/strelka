//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2018 Illumina, Inc.
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

/// \file
/// \author Chris Saunders
///


#include "AlleleReportInfoUtil.hh"
#include "blt_common/ref_context.hh"
#include "blt_util/seq_util.hh"
#include "starling_common/pos_basecall_buffer.hh"
#include "starling_common/readMappingAdjustmentUtil.hh"

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

static
void
set_repeat_info(
    const IndelKey& indelKey,
    const reference_contig_segment& ref,
    AlleleReportInfo& indelReportInfo)
{
    // Encoded lengths using fixed number of Zev-Lempel 1977 keywords. See http://www.lptmc.jussieu.fr/user/lesne/PRE-Short.pdf :
    indelReportInfo.contextCompressability = computeContextCompressability(ref, indelKey.pos, indelKey.right_pos(), 5);

    if (! ((indelReportInfo.it == SimplifiedIndelReportType::INSERT) ||
           (indelReportInfo.it == SimplifiedIndelReportType::DELETE) ||
           (indelReportInfo.it == SimplifiedIndelReportType::SWAP))) return;

    unsigned insert_repeat_count(0);
    unsigned delete_repeat_count(0);

    if      (indelReportInfo.it == SimplifiedIndelReportType::INSERT)
    {
        get_seq_repeat_unit(indelKey.insert_seq(), indelReportInfo.repeatUnit, insert_repeat_count);
    }
    else if (indelReportInfo.it == SimplifiedIndelReportType::DELETE)
    {
        std::string deletedSeq;
        copy_ref_subseq(ref,indelKey.pos, indelKey.pos+indelKey.delete_length(), deletedSeq);

        get_seq_repeat_unit(deletedSeq, indelReportInfo.repeatUnit, delete_repeat_count);
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

        indelReportInfo.repeatUnit=insert_ru;
    }
    else
    {
        assert(false && "Unexpected indel type");
    }
    indelReportInfo.repeatUnitLength=indelReportInfo.repeatUnit.size();

    // count repeats in contextual sequence:
    unsigned indel_context_repeat_count(0);
    {
        const pos_t indel_begin_pos(indelKey.pos);
        const pos_t indel_end_pos(indelKey.right_pos());
        const int repeat_unit_size(static_cast<int>(indelReportInfo.repeatUnit.size()));

        // count upstream repeats:
        for (pos_t i(indel_begin_pos-repeat_unit_size); i>=0; i-=repeat_unit_size)
        {
            bool is_repeat(true);
            for (int j(0); j<repeat_unit_size; ++j)
            {
                if (ref.get_base(i+j) != indelReportInfo.repeatUnit[j])
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
                if (ref.get_base(i+j) != indelReportInfo.repeatUnit[j])
                {
                    is_repeat = false;
                    break;
                }
            }
            if (! is_repeat) break;
            indel_context_repeat_count += 1;
        }
    }

    indelReportInfo.refRepeatCount = indel_context_repeat_count+delete_repeat_count;
    indelReportInfo.indelRepeatCount = indel_context_repeat_count+insert_repeat_count;
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
    indelReportInfo.interruptedHomopolymerLength = getInterruptedHomopolymerLength(indel_begin_pos - 1, ref);
    indelReportInfo.interruptedHomopolymerLength = std::max(indelReportInfo.interruptedHomopolymerLength,
                                                            getInterruptedHomopolymerLength(indel_begin_pos, ref));
    if (indel_begin_pos != indel_end_pos)
    {
        indelReportInfo.interruptedHomopolymerLength = std::max(indelReportInfo.interruptedHomopolymerLength,
                                                                getInterruptedHomopolymerLength(indel_end_pos - 1, ref));
        indelReportInfo.interruptedHomopolymerLength = std::max(indelReportInfo.interruptedHomopolymerLength,
                                                                getInterruptedHomopolymerLength(indel_end_pos, ref));
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

    // This term formally had a prior on incorrect mapping (meaning the prior of incorrectly mapping a read when the
    // read is randomly positioned on the genome). This term is effectively 1 so it is approximated out in the current
    // version.
    ReadPathScores::score_t pprob_readIncorrectlyMapped =
        getIncorrectMappingLogLikelihood(dopt, is_tier2_pass, path_lnp.nonAmbiguousBasesInRead);
    pprob.ref     = path_lnp.ref     + dopt.correctMappingLogPrior + allele_lnprior;
    pprob.indel   = path_lnp.indel   + dopt.correctMappingLogPrior + allele_lnprior;

    if (is_use_alt_indel)
    {
        for (const auto& val : path_lnp.alt_indel)
        {
            pprob.alt_indel.push_back(std::make_pair(val.first,(val.second + dopt.correctMappingLogPrior + allele_lnprior)));
        }
    }

    double scale(std::max(pprob_readIncorrectlyMapped,std::max(pprob.ref,pprob.indel)));

    if (is_use_alt_indel)
    {
        for (const auto& val : pprob.alt_indel)
        {
            if (scale < val.second) scale = val.second;
        }
    }

    pprob_readIncorrectlyMapped = std::exp(pprob_readIncorrectlyMapped-scale);
    pprob.ref = std::exp(pprob.ref-scale);
    pprob.indel = std::exp(pprob.indel-scale);
    if (is_use_alt_indel)
    {
        for (auto& val : pprob.alt_indel)
        {
            val.second = std::exp((val.second)-scale);
        }
    }


    double sum(pprob_readIncorrectlyMapped+pprob.ref+pprob.indel);
    if (is_use_alt_indel)
    {
        for (const auto& val : pprob.alt_indel)
        {
            sum += val.second;
        }
    }

    // pprob_readIncorrectlyMapped /= sum; /// no point in normalizing this if we aren't adding it back into pprob
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
    const bool isUseTier2Data,
    const bool isUseAltIndel,
    AlleleSampleReportInfo& isri)
{
    // get read info:
    {
        unsigned n_subscore_reads(0);

        for (const auto& val : indelSampleData.read_path_lnp)
        {
            const ReadPathScores& path_lnp(val.second);

            // optionally skip tier2 data:
            if ((! isUseTier2Data) && (! path_lnp.is_tier1_read)) continue;

            const ReadPathScores pprob(indel_lnp_to_pprob(dopt, path_lnp, isUseTier2Data, isUseAltIndel));

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

        if (isUseTier2Data)
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
        isri.indelLocusDepth=spi.calls.size();
        isri.mapqTracker=spi.mapqTracker;
    }
}

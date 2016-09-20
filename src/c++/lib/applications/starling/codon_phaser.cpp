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
 *  Created on: Sep 10, 2013
 *  Author: Morten Kallberg
 */

#include "codon_phaser.hh"

#include <array>
#include <functional>
#include <vector>

//#define DEBUG_CODON



#ifdef DEBUG_CODON
#include "blt_util/log.hh"
#endif


void
Codon_phaser::
process(std::unique_ptr<GermlineSiteLocusInfo> locusPtr)
{
    if (_opt.isUseCodonPhaser)
    {
        // make one pass through to determine if this site will be buffered
        bool isBufferSite(false);
        const unsigned sampleCount(_sampleBlocks.size());
        for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
        {
            auto& block(_sampleBlocks[sampleIndex]);
            if (is_phasable_locus(*locusPtr, sampleIndex))
            {
                // extending block with het call, update block_end position
                isBufferSite=true;
                if (! block.is_in_block())
                {
                    block.start = locusPtr->pos;
                }
                block.end = locusPtr->pos;
                block.het_count++;

            }
            else if (block.is_in_block())
            {
                if ((locusPtr->pos - block.end + 1) < _opt.phasing_window)
                {
                    // extending block with non-het call based on the phasing range
                    isBufferSite = true;
                }
                else
                {
                    // past phasing window, no phasing oppurtunity in this sample:
                    block.clear();
                }
            }
        }

        if (isBufferSite)
        {
            _buffer.push_back(std::move(locusPtr));
            for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
            {
                auto& block(_sampleBlocks[sampleIndex]);
                if (block.het_count == 2)
                {
                    phaseRecords(sampleIndex);
                    block.start = block.end;
                    block.het_count = 1;
                }
            }
            return;
        }
        else
        {
            output_buffer();
        }
    }
    _sink->process(std::move(locusPtr));
}



void
Codon_phaser::
process(std::unique_ptr<GermlineIndelLocusInfo> locusPtr)
{
    // the Codon phaser can't work across indels, so flush any in-progress phasing
    if (_opt.isUseCodonPhaser)
    {
        output_buffer();
    }
    _sink->process(std::move(locusPtr));
}



void
Codon_phaser::
flush_impl()
{
    if (_opt.isUseCodonPhaser)
    {
        output_buffer();
    }
}


void
Codon_phaser::
output_buffer()
{
    if (not isBuffer()) return;

    for (auto& siteLocus : _buffer)
    {
        _sink->process(std::move(siteLocus));
    }

    for (auto& sampleBlock : _sampleBlocks)
    {
        sampleBlock.clear();
    }
    _buffer.clear();
}



bool
Codon_phaser::
evaluatePhasingConsistency(
    const BlockReadObservations& blockObs,
    const std::array<std::pair<std::string,allele_observations>, 2>& max_alleles)
{
    // some ad hoc metrics to measure consistency with diploid model
    const int allele_sum = max_alleles[0].second.count() + max_alleles[1].second.count();
    const float max_allele_frac = (static_cast<float>(allele_sum))/blockObs.total_reads;
    const float relative_allele_frac = static_cast<float>(max_alleles[1].second.count())/max_alleles[0].second.count();

#ifdef DEBUG_CODON
    log_os << "max alleles sum " << allele_sum << "\n";
    log_os << "max alleles frac " << max_allele_frac << "\n";
    log_os << "relative_allele_frac " << relative_allele_frac << "\n";
#endif

    if (max_allele_frac<0.8)
    {
#ifdef DEBUG_CODON
        log_os << __FUNCTION__ << "; non-diploid\n";
#endif
        // non-diploid?
        return false;
    }

    if (relative_allele_frac<0.5)
    {
        // allele imbalance?
#ifdef DEBUG_CODON
        log_os << __FUNCTION__ << "; allele imbalance\n";
#endif
        return false;
    }

    // sanity check that we have one het on each end of the block:
    if (max_alleles[0].second.count() == 0) return false;
    if (max_alleles[1].second.count() == 0) return false;

    assert(! max_alleles[0].first.empty());
    assert(max_alleles[0].first.size() == max_alleles[1].first.size());
    if (max_alleles[0].first.front() == max_alleles[1].first.front()) return false;
    if (max_alleles[0].first.back() == max_alleles[1].first.back()) return false;

    return true;
}


struct PhaseEditInfo
{
    unsigned index;
    bool isFlip;
};


void
Codon_phaser::
create_phased_record(
    const unsigned sampleIndex,
    const BlockReadObservations& blockObs)
{
    const auto& block(_sampleBlocks[sampleIndex]);

    // sanity check do we have all record we expect or are there un-accounted no-calls
    if (block.get_length() > _buffer.size())
        return;

    if (blockObs.total_reads<10)
    {
        // some initial minimum conditions, look for at least 10 spanning reads support
        return;
    }

    assert(_buffer.size()>0);

    /// TODO STREL-125 TMP:
    assert(_buffer[0]->getSampleCount() == _sampleBlocks.size());

    //Decide if we accept the novel alleles, very hacky criteria for now
    //at least 80% of the reads have to support a diploid model
    //TODO unphased genotype corresponds to as many phased genotypes as there are permutations of
    //the observed alleles. Thus, for a given unphased genotyping G1, . . . ,Gn,
    //we need to to calculate probability of sampling a particular unphased genotype combinations
    //given a set of allele frequencies...
    // TODO: PASSing SNPs are getting filtered during phasing for allele balance issues. Looks like it would benefit from a
    // model-based approach.
    // < chr8  70364286    .   G   A   32  PASS    SNVSB=-2.5;SNVHPOL=4    GT:GQ:GQX:DP:DPF:AD 0/1:65:21:15:5:12,3
    // < chr8  70364287    .   A   G   11  PASS    SNVSB=-2.5;SNVHPOL=2    GT:GQ:GQX:DP:DPF:AD 0/1:44:18:16:3:14,2
    //---
    //> chr8  70364286    .   G   A   32  PhasingConflict SNVSB=-2.5;SNVHPOL=4    GT:GQ:GQX:DP:DPF:AD 0/1:65:21:15:5:12,3
    //> chr8  70364287    .   A   G   11  PhasingConflict SNVSB=-2.5;SNVHPOL=2    GT:GQ:GQX:DP:DPF:AD 0/1:44:18:16:3:14,2
    //static constexpr unsigned alleleCount(2);
    std::array<std::pair<std::string,allele_observations>, 2> max_alleles;
    for (const auto& obs : blockObs.observations)
    {
#ifdef DEBUG_CODON
        log_os << "obs:" << obs.first << "(" << obs.second.count() << ")\n";
#endif

        if (obs.second.count()>max_alleles[0].second.count())
        {
            max_alleles[1] = max_alleles[0];
            max_alleles[0] = obs;
        }
        if (obs.second.count()>max_alleles[1].second.count() && max_alleles[0].first!=obs.first)
        {
            max_alleles[1] = obs;
        }
    }

#ifdef DEBUG_CODON
    log_os << "max_1 " << max_alleles[0].first << "=" << max_alleles[0].second.count() << "\n";
    log_os << "max_2 " << max_alleles[1].first << "=" << max_alleles[1].second.count() << "\n";
    log_os << "buffer size " << _buffer.size() << "\n";
    log_os << "block length " << block.get_length() << "\n";
#endif

    // where does this block start in the buffer:
    const int blockStartOffset(block.start - bufferStartPos());
    assert(blockStartOffset >= 0);
    const unsigned blockLength(block.get_length());

    bool isPhasingConsistent(evaluatePhasingConsistency(blockObs, max_alleles));

    // edit info:
    std::vector<PhaseEditInfo> phaseEditInfo;
    pos_t phaseSetId(-1);

    if (isPhasingConsistent)
    {
        const std::string& phasedHap0(max_alleles[0].first);
        const std::string& phasedHap1(max_alleles[1].first);

        bool isMatchedAtBlockStart(true);

        for (unsigned blockIndex(0); blockIndex < blockLength; blockIndex++)
        {
            const GermlineSiteLocusInfo& siteLocus(*_buffer.at(blockStartOffset + blockIndex));
            if (!is_phasable_locus(siteLocus, sampleIndex)) continue;

            const auto& sampleInfo(siteLocus.getSample(sampleIndex));
            const auto& maxGt(sampleInfo.max_gt());

            const uint8_t allele0Index(maxGt.getAllele0Index());
            const uint8_t allele1Index(maxGt.getAllele1Index());

            const auto& siteAlleles(siteLocus.getSiteAlleles());

            auto alleleIndexToBaseIndex = [&](const uint8_t alleleIndex)
            {
                if (alleleIndex == 0) return siteLocus.refBaseIndex;
                return static_cast<uint8_t>(siteAlleles[alleleIndex - 1].baseIndex);
            };

            const uint8_t base0Index(alleleIndexToBaseIndex(allele0Index));
            const uint8_t base1Index(alleleIndexToBaseIndex(allele1Index));

            const uint8_t phasedBase0Index(base_to_id(phasedHap0[blockIndex]));
            const uint8_t phasedBase1Index(base_to_id(phasedHap1[blockIndex]));

            bool isMatched(true);
            if ((base0Index == phasedBase0Index) and (base1Index == phasedBase1Index))
            {
                isMatched = true;
            }
            else if ((base0Index == phasedBase1Index) and (base1Index == phasedBase0Index))
            {
                isMatched = false;
            }
            else
            {
                isPhasingConsistent = false;
                break;
            }

            bool isFlip(false);
            if (phaseSetId == -1)
            {
                if (sampleInfo.phaseSetId != -1)
                {
                    phaseSetId = sampleInfo.phaseSetId;
                }
                else
                {
                    phaseSetId = (siteLocus.pos + 1);
                }

                isMatchedAtBlockStart = isMatched;
            }
            else
            {
                isFlip = (isMatched != isMatchedAtBlockStart);
            }
            phaseEditInfo.push_back(PhaseEditInfo({blockIndex, isFlip}));
        }
    }

    if (isPhasingConsistent)
    {
        // if everything went wellgo back through one more time and make all edits:
        for (const auto& editInfo : phaseEditInfo)
        {
            LocusInfo& locus(*_buffer.at(blockStartOffset + editInfo.index));
            auto& sampleInfo(locus.getSample(sampleIndex));
            sampleInfo.phaseSetId = phaseSetId;
            sampleInfo.max_gt().setPhased(editInfo.isFlip);
        }
    }
    else
    {
        // should a phasing conflict filter be provided, or should be just leave the
        // variants unphased/unfiltered? Choosing the latter option here.
#if 0
        for (unsigned blockIndex(0); blockIndex < blockLength; blockIndex++)
        {
            GermlineSiteLocusInfo& siteLocus(*_buffer.at(blockStartOffset + blockIndex));
            if (! is_phasable_site(siteLocus, sampleIndex)) continue;
            auto& sampleInfo(siteLocus.getSample(sampleIndex));
            sampleInfo.filters.set(GERMLINE_VARIANT_VCF_FILTERS::PhasingConflict);
        }
#endif
    }
}



// makes the phased VCF record from the buffered sites list
void
Codon_phaser::
phaseRecords(const unsigned sampleIndex)
{
    BlockReadObservations blockObs;
    collect_pileup_evidence(sampleIndex, blockObs);
    create_phased_record(sampleIndex, blockObs);
}



void
Codon_phaser::
collect_pileup_evidence(
    const unsigned sampleIndex,
    BlockReadObservations& blockObs)
{
    const auto& block(_sampleBlocks[sampleIndex]);
    const pos_basecall_buffer& bc_buff(_basecallBuffers[sampleIndex]);

    // build quick pileup index over phase range:
    std::vector<const snp_pos_info*> spi;
    for (int blockPos(block.start); blockPos<=block.end; ++blockPos)
    {
        const snp_pos_info& pi(bc_buff.get_pos(blockPos));
        spi.push_back(&pi);
    }

    const unsigned blockWidth(spi.size());
    std::vector<int> callOffset(blockWidth,0);

    /// traces individual read fragments from the pileup structure within the phasing block range
    ///
    /// Low detail summary:
    /// If translating aligned reads into pileup columns is thought of as a sort of matrix
    /// transpose, this function is trying to invert the transposition back to a (partial)
    /// read. A naive pileup structure would not support this, but starling pileup information
    /// has been supplemented with a few extra bits that allows this reconstruction.
    ///
    /// Given a start offset within the phasing block (startBlockIndex), return a tuple composed of:
    ///     1. A bool indicating whether the read is good (ie. complete to the end of the
    ///        phasing block and all basecalls pass filter) This return val is really only
    ///        useful for startBlockIndex=0
    ///     2. The reconstructed read fragment
    ///
    /// External calls should always provide startBlockIndex=0, other values of this argument
    /// are used by internal recursive calls to the function. If the return bool indicates a
    /// good read, this function effectively returns a single allele count for the phaser
    ///
    /// Important implementation detail: this function mutates the external vector 'callOffset'
    /// over successive calls to progressively jump to the offset of the next read
    ///
    /// isFirstBaseCallFromMatchSeg is set for any basecall starting a continuous matching sequence in one read
    /// isLastBaseCallFromMatchSeg is set for any basecall ending a continuous matching sequence in one read
    // can't use auto for the return value here b/c of recursion:
    std::function<std::pair<bool,std::string>(const unsigned)> tracePartialRead =
        [&](const unsigned startBlockIndex)
    {
        bool isPass(true);
        std::string readFrag;
        for (unsigned blockIndex(startBlockIndex); blockIndex<blockWidth; ++blockIndex)
        {
            const snp_pos_info& pi(*spi[blockIndex]);
            while (true)
            {
                const base_call& bc(pi.calls.at(callOffset[blockIndex]));
                if ((! bc.isFirstBaseCallFromMatchSeg) || (blockIndex == startBlockIndex)) break;
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wuninitialized"
                tracePartialRead(blockIndex);
#pragma clang diagnostic pop
            }
            const base_call& bc(pi.calls.at(callOffset[blockIndex]));
            // this represents the 'ordinary' advance of callOffset for a non-interrupted read segment:
            callOffset[blockIndex]++;
            if (bc.isLastBaseCallFromMatchSeg && ((blockIndex+1) < blockWidth))
            {
                isPass=false;
                break;
            }
            if (bc.is_call_filter) isPass=false;
            readFrag += id_to_base(bc.base_id);
        }
        return std::make_pair(isPass,readFrag);
    };

    // analyze as virtual reads -- to do so treat the first pileup column as a privileged reference point:
    const snp_pos_info& beginPi(*spi[0]);
    const unsigned n_calls(beginPi.calls.size());
    for (unsigned callIndex(0); callIndex<n_calls; ++callIndex)
    {
        const bool is_fwd_strand(beginPi.calls[callIndex].is_fwd_strand);
        std::pair<bool,std::string> result = tracePartialRead(0);
#ifdef DEBUG_CODON
        log_os << __FUNCTION__ << ": callOffset:";
        for (const auto off : callOffset)
        {
            log_os << "\t" << off;
        }
        log_os << "\n";
        log_os << __FUNCTION__ << ": callIndex, isPass, substr: " << callIndex << " " << result.first << " " << result.second << "\n";
#endif

        if (! result.first) continue;

        if (is_fwd_strand)
        {
            blockObs.observations[result.second].fwd++;
        }
        else
        {
            blockObs.observations[result.second].rev++;
        }
        blockObs.total_reads++;
    }

    blockObs.total_reads_unused=n_calls-blockObs.total_reads;
}

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
        if (locusPtr->getActiveRegionId() != _activeRegionId)
        {
            
        }


//        LocusInfo& locus(*_buffer.at(blockStartOffset + editInfo.index));
//        auto& sampleInfo(locus.getSample(sampleIndex));
//        sampleInfo.phaseSetId = phaseSetId;
//        sampleInfo.max_gt().setPhased(editInfo.isFlip);
    }
    _sink->process(std::move(locusPtr));
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

    _buffer.clear();
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
create_phased_record(
        const unsigned sampleIndex)
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


void
Codon_phaser::
flush_impl()
{
    if (_opt.isUseCodonPhaser)
    {
        output_buffer();
    }
}
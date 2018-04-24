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
/// \author Sangtae Kim
///

#include "boost/test/unit_test.hpp"

#include "starling_common/ActiveRegionDetector.hh"
#include "starling_common/CandidateSnvBuffer.hh"
#include "test/testIndelBuffer.hh"


/// Test whether ActiveRegionDetector correctly handles ploidy
static int ploidyTest(const unsigned ploidy)
{
    assert(ploidy <= 4);
    reference_contig_segment ref;
    ref.seq() = "GATCTGT";
    const unsigned maxIndelSize(49);
    const unsigned sampleCount(1);
    const unsigned sampleIndex(0);
    const unsigned depth(50);

    TestIndelBuffer testBuffer(ref);
    CandidateSnvBuffer testSnvBuffer(sampleCount);

    ActiveRegionDetector activeRegionDetector(ref, testBuffer.getIndelBuffer(),
                                              testSnvBuffer, maxIndelSize, sampleCount, false);

    const auto snvPos = std::set<pos_t>({2, 4});

    pos_t refLength = (pos_t)ref.seq().length();

    // create 4 haplotypes with differing bases at positions 2 and 4
    // hap0 (no SNV): 20 reads
    // hap1 (SNV at 2): 13 haplotypes
    // hap2 (SNV at 4): 12 haplotypes
    // hap3 (SNV at 6): 5 haplotypes

    // Given ploidy n, n haplotypes will be selected,
    // and thus #candidate SNVs will be n-1
    for (unsigned alignId(0); alignId < depth; ++alignId)
    {
        bool isForwardStrand = ((alignId % 2) == 0);
        activeRegionDetector.getReadBuffer(sampleIndex).setAlignInfo(
            alignId, sampleIndex, INDEL_ALIGN_TYPE::GENOME_TIER1_READ, isForwardStrand);

        bool isSnvAtPos2(false);
        bool isSnvAtPos4(false);
        bool isSnvAtPos6(false);
        if (alignId < 20)
        {
            // no SNV
        }
        else if (alignId < 33)
        {
            isSnvAtPos2 = true;
        }
        else if (alignId < 45)
        {
            isSnvAtPos4 = true;
        }
        else
        {
            isSnvAtPos6 = true;
        }

        for (pos_t pos(0); pos<refLength; ++pos)
        {

            if ((isSnvAtPos2 && (pos == 2))
                || (isSnvAtPos4 && (pos == 4))
                || (isSnvAtPos6 && (pos == 6)))
            {
                // SNV position
                activeRegionDetector.getReadBuffer(sampleIndex).insertMismatch(alignId, pos, 'A');
            }
            else
            {
                // No SNV
                activeRegionDetector.getReadBuffer(sampleIndex).insertMatch(alignId, pos);
            }
        }
    }

    // Create and process active regions
    for (pos_t pos(0); pos<refLength; ++pos)
    {
        activeRegionDetector.updateSamplePloidy(sampleIndex, pos, ploidy);
        activeRegionDetector.updateEndPosition(pos);
    }
    activeRegionDetector.clear();

    // count #candidate SNVs
    unsigned numCandidateSNVs(0);
    for (pos_t pos(0); pos<refLength; ++pos)
    {
        if (testSnvBuffer.isCandidateSnvAnySample(pos, 'A'))
            ++numCandidateSNVs;
    }

    return numCandidateSNVs;
}

BOOST_AUTO_TEST_SUITE( test_activeRegionPloidy )

BOOST_AUTO_TEST_CASE( test_activeRegionPloidy1To3 )
{
    for (unsigned ploidy(1); ploidy<4; ++ploidy)
    {
        BOOST_REQUIRE_EQUAL(ploidyTest(ploidy), ploidy-1);
    }
}

BOOST_AUTO_TEST_SUITE_END()

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

const unsigned maxIndelSize = 49;


BOOST_AUTO_TEST_SUITE( test_activeRegion )

// checks whether positions with consistent mismatches are marked as polymorphic sites
BOOST_AUTO_TEST_CASE( test_multiSampleMMDF )
{
    reference_contig_segment ref;
    ref.seq() = "GATCTGT";
    const int sampleCount = 3;
    const int depth = 50;

    TestIndelBuffer testBuffer(ref);
    CandidateSnvBuffer testSnvBuffer(sampleCount);

    ActiveRegionDetector activeRegionDetector(ref, testBuffer.getIndelBuffer(),
                                              testSnvBuffer, maxIndelSize, sampleCount, false);

    const auto snvPos = std::set<pos_t>({2, 4, 5});

    pos_t refLength = (pos_t)ref.seq().length();

    // fake reading reads
    for (int alignId=0; alignId < depth*sampleCount; ++alignId)
    {
        unsigned sampleIndex = alignId % sampleCount;
        bool isForwardStrand = ((alignId % 4) == 0) or ((alignId % 4) == 3);
        activeRegionDetector.getReadBuffer(sampleIndex).setAlignInfo(alignId, sampleIndex, INDEL_ALIGN_TYPE::GENOME_TIER1_READ, isForwardStrand);
        for (pos_t pos(0); pos<refLength; ++pos)
        {
            // only sample 0 has mismatches
            // alternative allele frequency 0.5
            if (sampleIndex != 0 or ((alignId/sampleCount) % 2)
                or (snvPos.find(pos) == snvPos.end()))
            {
                // No SNV
                activeRegionDetector.getReadBuffer(sampleIndex).insertMatch(alignId, pos);
            }
            else
            {
                // SNV position
                activeRegionDetector.getReadBuffer(sampleIndex).insertMismatch(alignId, pos, 'A');
            }
        }
    }

    for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
    {
        for (pos_t pos(0); pos<refLength; ++pos)
        {
            activeRegionDetector.updateEndPosition(pos);
        }
        activeRegionDetector.clear();
    }

    // check if isCandidateSnv are correctly set
    for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
    {
        for (pos_t pos(0); pos<refLength; ++pos)
        {
            if (snvPos.find(pos) != snvPos.end())
            {
                // SNV
                BOOST_REQUIRE_EQUAL(testSnvBuffer.isCandidateSnvAnySample(pos, 'A'), true);
            }
            else
            {
                // No SNV
                BOOST_REQUIRE_EQUAL(testSnvBuffer.isCandidateSnvAnySample(pos, 'A'), false);
            }
        }
    }
}

// Checks whether an indel is correctly confirmed in active regions
BOOST_AUTO_TEST_CASE( test_indelCandidacy )
{
    reference_contig_segment ref;
    ref.seq() = "TCTCT";

    const int sampleCount = 1;
    const unsigned sampleIndex = 0;

    TestIndelBuffer testBuffer(ref);
    CandidateSnvBuffer testSnvBuffer(sampleCount);

    ActiveRegionDetector detector(ref, testBuffer.getIndelBuffer(),
                                  testSnvBuffer, maxIndelSize, sampleCount, false);

    const int depth = 50;

    const pos_t indelPos = 2;
    auto indelKey = IndelKey(indelPos, INDEL::INDEL, 0, "AG");

    pos_t refLength = (pos_t)ref.seq().length();

    // fake reading reads
    for (int alignId=0; alignId < depth; ++alignId)
    {
        bool isForwardStrand = ((alignId % 4) == 0) or ((alignId % 4) == 3);
        detector.getReadBuffer(sampleIndex).setAlignInfo(alignId, sampleIndex, INDEL_ALIGN_TYPE::GENOME_TIER1_READ, isForwardStrand);
        for (pos_t pos(0); pos<refLength; ++pos)
        {
            detector.getReadBuffer(sampleIndex).insertMatch(alignId, pos);

            if (pos == indelPos && (alignId % 2))
            {
                IndelObservation indelObservation;

                IndelObservationData indelObservationData;
                indelObservation.key = indelKey;
                indelObservationData.id = alignId;
                indelObservationData.iat = INDEL_ALIGN_TYPE::GENOME_TIER1_READ;
                indelObservation.data = indelObservationData;

                detector.getReadBuffer(sampleIndex).insertIndel(indelObservation);
            }
        }
    }

    for (pos_t pos(0); pos<refLength; ++pos)
    {
        detector.updateEndPosition(pos);
    }
    detector.clear();

    const auto itr(testBuffer.getIndelBuffer().getIndelIter(indelKey));
    BOOST_REQUIRE_EQUAL(itr->second.isDiscoveredInActiveRegion(), true);
}


BOOST_AUTO_TEST_CASE( test_jumpingPositions )
{
    const pos_t startPositions[] = {100, 7000};
    const pos_t snvOffsets[] = {40, 42};

    std::string refSeq(10000, 'A');
    for (auto startPosition : startPositions)
    {
        for (auto snvOffset : snvOffsets)
        {
            // to avoid hpol filter
            refSeq[startPosition+snvOffset - 1] = 'T';
            refSeq[startPosition+snvOffset + 1] = 'T';
        }
    }

    reference_contig_segment ref;
    ref.seq() = refSeq;

    const int sampleCount = 1;
    const unsigned sampleIndex = 0;

    TestIndelBuffer testBuffer(ref);
    CandidateSnvBuffer testSnvBuffer(sampleCount);

    ActiveRegionDetector detector(ref, testBuffer.getIndelBuffer(),
                                  testSnvBuffer, maxIndelSize, sampleCount, false);

    // fake reading reads
    const int depth = 50;
    const unsigned readLength = 100;


    for (auto startPosition : startPositions)
    {
        pos_t endPosition = startPosition + readLength;
        for (int alignId=0; alignId < depth; ++alignId)
        {
            bool isForwardStrand = ((alignId % 4) == 0) or ((alignId % 4) == 3);
            detector.getReadBuffer(sampleIndex).setAlignInfo(alignId, sampleIndex, INDEL_ALIGN_TYPE::GENOME_TIER1_READ, isForwardStrand);
            for (pos_t pos(startPosition); pos<endPosition; ++pos)
            {
                // SNVs at startPosition+10 and startPosition+12
                auto isSnvPosition = (pos == startPosition+snvOffsets[0]) || (pos == startPosition+snvOffsets[1]);
                if ((alignId % 2) && isSnvPosition)
                {
                    detector.getReadBuffer(sampleIndex).insertMismatch(alignId, pos, 'G');
                }
                else
                {
                    detector.getReadBuffer(sampleIndex).insertMatch(alignId, pos);
                }
            }
        }

        for (pos_t pos(startPosition); pos<endPosition; ++pos)
        {
            detector.updateEndPosition(pos);
        }
        detector.clear();

        // check if isCandidateSnv are correctly set
        BOOST_REQUIRE_EQUAL(testSnvBuffer.isCandidateSnvAnySample(startPosition + snvOffsets[0], 'G'), true);
        BOOST_REQUIRE_EQUAL(testSnvBuffer.isCandidateSnvAnySample(startPosition + snvOffsets[1], 'G'), true);
    }
}

// Checks whether an indel is left-shifted
BOOST_AUTO_TEST_CASE( test_leftShiftIndel )
{
    reference_contig_segment ref;
    ref.seq() = "GTCC";

    const unsigned sampleCount = 1;
    const unsigned sampleIndex = 0;

    TestIndelBuffer testBuffer(ref);
    CandidateSnvBuffer testSnvBuffer(sampleCount);

    ActiveRegionDetector detector(ref, testBuffer.getIndelBuffer(),
                                  testSnvBuffer, maxIndelSize, sampleCount, false);

    const int depth = 50;

    const pos_t indelPos = 2;
    auto indelKey = IndelKey(indelPos, INDEL::INDEL, 0, "ATAT");

    pos_t refLength = (pos_t)ref.seq().length();

    // fake reading reads
    for (int alignId=0; alignId < depth; ++alignId)
    {
        bool isForwardStrand = ((alignId % 4) == 0) or ((alignId % 4) == 3);
        detector.getReadBuffer(sampleIndex).setAlignInfo(alignId, sampleIndex, INDEL_ALIGN_TYPE::GENOME_TIER1_READ, isForwardStrand);
        for (pos_t pos(0); pos<refLength; ++pos)
        {
            detector.getReadBuffer(sampleIndex).insertMatch(alignId, pos);

            if (pos == indelPos && (alignId % 2))
            {
                IndelObservation indelObservation;

                IndelObservationData indelObservationData;
                indelObservation.key = indelKey;
                indelObservationData.id = alignId;
                indelObservationData.iat = INDEL_ALIGN_TYPE::GENOME_TIER1_READ;
                indelObservation.data = indelObservationData;

                detector.getReadBuffer(sampleIndex).insertIndel(indelObservation);
            }
        }
    }

    for (pos_t pos(0); pos<refLength; ++pos)
    {
        detector.updateEndPosition(pos);
    }
    detector.clear();

    // check if the indel is shifted 1 base to the left
    auto leftShiftedIndelKey = IndelKey(indelPos-1, INDEL::INDEL, 0, "TATA");

    const auto itr(testBuffer.getIndelBuffer().getIndelIter(leftShiftedIndelKey));
    BOOST_REQUIRE_EQUAL(itr->second.isDiscoveredInActiveRegion(), true);
}


BOOST_AUTO_TEST_SUITE_END()

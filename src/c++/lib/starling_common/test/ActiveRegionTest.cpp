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

///
/// \author Sangtae Kim
///

#include <starling_common/IndelBuffer.hh>
#include <starling_common/ActiveRegionDetector.hh>
#include "boost/test/unit_test.hpp"

BOOST_AUTO_TEST_SUITE( test_activeRegion )

typedef std::unique_ptr<IndelBuffer> IndelBufferPtr;

struct TestIndelBuffer
{
    explicit
    TestIndelBuffer(
        const reference_contig_segment& ref)
    {
        // fake starling options
        _opt.is_user_genome_size = true;
        _opt.user_genome_size = ref.seq().size();

        const double maxDepth = 100.0;
        _doptPtr.reset(new starling_base_deriv_options(_opt,ref));

        _IndelBufferPtr.reset(new IndelBuffer(_opt, *_doptPtr, ref));
        _IndelBufferPtr->registerSample(depth_buffer(), depth_buffer(), maxDepth);
        _IndelBufferPtr->finalizeSamples();

    }

    IndelBuffer&
    getIndelBuffer()
    {
        return *_IndelBufferPtr;
    }

private:
    starling_base_options _opt;
    std::unique_ptr<starling_base_deriv_options> _doptPtr;
    std::unique_ptr<IndelBuffer> _IndelBufferPtr;
};


// checks whether positions with consistent mismatches are marked as polymorphic sites
BOOST_AUTO_TEST_CASE( test_multiSampleMMDF )
{
    reference_contig_segment ref;
    ref.seq() = "TCTGT";
    TestIndelBuffer testBuffer(ref);

    const unsigned maxIndelSize = 50;
    const int sampleCount = 3;
    const int depth = 50;

    ActiveRegionDetector detector(ref, testBuffer.getIndelBuffer(), maxIndelSize, sampleCount);

    const auto snvPos = std::set<pos_t>({0, 2, 3});

    pos_t refLength = (pos_t)ref.seq().length();

    // fake reading reads
    for (int alignId=0; alignId < depth*sampleCount; ++alignId)
    {
        unsigned sampleId = alignId % sampleCount;
        detector.setAlignInfo(alignId, sampleId, INDEL_ALIGN_TYPE::GENOME_TIER1_READ);
        for (pos_t pos(0); pos<refLength-1; ++pos)
        {
            // only sample 1 has mismatches
            // alternative allele frequency 0.5
            if (sampleId != 1 or ((alignId/sampleCount) % 2)
                or (snvPos.find(pos) == snvPos.end()))
            {
                // No SNV
                detector.insertMatch(alignId, pos);
            }
            else
            {
                // SNV position
                detector.insertMismatch(alignId, pos, 'A');
            }
        }
    }

    for (pos_t pos(0); pos<refLength-1; ++pos)
    {
        detector.updateEndPosition(pos, false);
    }
    detector.updateEndPosition(refLength-1, true);


    // check if polySites are correctly set

    for (unsigned sampleId(0); sampleId<sampleCount; ++sampleId)
    {
        for (pos_t pos(0); pos<refLength-1; ++pos)
        {
            if ((sampleId == 1) and (snvPos.find(pos) != snvPos.end()))
            {
                // SNV
                BOOST_REQUIRE_EQUAL(detector.isPolymorphicSite(sampleId, pos), true);
            }
            else
            {
                // No SNV
                BOOST_REQUIRE_EQUAL(detector.isPolymorphicSite(sampleId, pos), false);
            }
        }
    }
}

// Checks whether an indel is correctly confirmed in active regions
BOOST_AUTO_TEST_CASE( test_indelCandidacy )
{
    reference_contig_segment ref;
    ref.seq() = "TCTCT";

    TestIndelBuffer testBuffer(ref);

    const unsigned maxIndelSize = 50;
    const unsigned sampleCount = 1;
    ActiveRegionDetector detector(ref, testBuffer.getIndelBuffer(), maxIndelSize, sampleCount);

    const int sampleId = 0;
    const int depth = 50;

    const pos_t indelPos = 2;
    auto indelKey = IndelKey(indelPos, INDEL::INDEL, 0, "AG");

    pos_t refLength = (pos_t)ref.seq().length();

    // fake reading reads
    for (int alignId=0; alignId < depth; ++alignId)
    {
        detector.setAlignInfo(alignId, sampleId, INDEL_ALIGN_TYPE::GENOME_TIER1_READ);
        for (pos_t pos(0); pos<refLength-1; ++pos)
        {
            detector.insertMatch(alignId, pos);

            if (pos == indelPos && (alignId % 2))
            {
                IndelObservation indelObservation;

                IndelObservationData indelObservationData;
                indelObservation.key = indelKey;
                indelObservationData.id = alignId;
                indelObservationData.iat = INDEL_ALIGN_TYPE::GENOME_TIER1_READ;
                indelObservation.data = indelObservationData;

                detector.insertIndel(sampleId, indelObservation);
            }
        }
    }

    for (pos_t pos(0); pos<refLength-1; ++pos)
    {
        detector.updateEndPosition(pos, false);
    }
    detector.updateEndPosition(refLength-1, true);

    const auto itr(testBuffer.getIndelBuffer().getIndelIter(indelKey));
    BOOST_REQUIRE_EQUAL(itr->second.isConfirmedInActiveRegion, true);
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

    TestIndelBuffer testBuffer(ref);

    const unsigned maxIndelSize = 50;
    const unsigned sampleCount = 1;
    ActiveRegionDetector detector(ref, testBuffer.getIndelBuffer(), maxIndelSize, sampleCount);

    // fake reading reads
    const int depth = 50;
    const int sampleId = 0;
    const unsigned readLength = 100;


    for (auto startPosition : startPositions)
    {
        pos_t endPosition = startPosition + readLength;
        for (int alignId=0; alignId < depth; ++alignId)
        {
            detector.setAlignInfo(alignId, sampleId, INDEL_ALIGN_TYPE::GENOME_TIER1_READ);
            for (pos_t pos(startPosition); pos<endPosition; ++pos)
            {
                // SNVs at startPosition+10 and startPosition+12
                auto isSnvPosition = (pos == startPosition+snvOffsets[0]) || (pos == startPosition+snvOffsets[1]);
                if ((alignId % 2) && isSnvPosition)
                {
                    detector.insertMismatch(alignId, pos, 'G');
                }
                else
                {
                    detector.insertMatch(alignId, pos);
                }
            }
        }

        for (pos_t pos(startPosition); pos<endPosition; ++pos)
        {
            detector.updateEndPosition(pos, false);
            detector.updateStartPosition(pos - (readLength + maxIndelSize));
        }

        // check if polySites are correctly set
        BOOST_REQUIRE_EQUAL(detector.isPolymorphicSite(sampleId, startPosition+snvOffsets[0]), true);
        BOOST_REQUIRE_EQUAL(detector.isPolymorphicSite(sampleId, startPosition+snvOffsets[1]), true);
    }
}


BOOST_AUTO_TEST_SUITE_END()

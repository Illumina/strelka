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
/// \author Sangtae Kim
///

#include <blt_util/reference_contig_segment.hh>
#include <starling_common/ReferenceRepeatFinder.hh>
#include "boost/test/unit_test.hpp"

BOOST_AUTO_TEST_SUITE( test_referenceRepeatFinder )

// Checks whether anchors are selected outside of repeat regions
BOOST_AUTO_TEST_CASE( test_findingAnchors )
{

    reference_contig_segment ref;

    // all the bases are repeat except AATG
    ref.seq() = "TATATACCCCCAATGAAAAA";
    const pos_t length(ref.seq().length());

    const unsigned maxRepeatLength(10u);
    const unsigned maxBufferSize(length);

    // AA is not considered as repeat because the size is less than minRepeatSpan
    const pos_t minRepeatSpan(3u);

    ReferenceRepeatFinder repeatFinder(ref, maxRepeatLength, maxBufferSize, minRepeatSpan);

    repeatFinder.initRepeatSpan(0u);

    for (pos_t pos(1u); pos<length; ++pos)
    {
        repeatFinder.updateRepeatSpan(pos);
    }

    for (pos_t pos(0); pos<length; ++pos)
    {
        if (pos >= 11 and pos <= 14)
            BOOST_CHECK(repeatFinder.isAnchor(pos));
        else
            BOOST_CHECK(not repeatFinder.isAnchor(pos));
    }
}

BOOST_AUTO_TEST_SUITE_END()

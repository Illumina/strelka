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

#include "boost/test/unit_test.hpp"

#include "blt_util/IntegerLogCompressor.hh"



BOOST_AUTO_TEST_SUITE( test_integerLogCompressor )

BOOST_AUTO_TEST_CASE( test_compression )
{
    static const unsigned bitCount(3);
    for (unsigned i(0); i<8; ++i)
    {
        BOOST_REQUIRE_EQUAL(compressInt(i,bitCount),i);
    }

    for (unsigned i(8); i<10; ++i)
    {
        BOOST_REQUIRE_EQUAL(compressInt(i,bitCount),9u);
    }
    for (unsigned i(10); i<12; ++i)
    {
        BOOST_REQUIRE_EQUAL(compressInt(i,bitCount),10u);
    }

    for (unsigned i(16); i<20; ++i)
    {
        BOOST_REQUIRE_EQUAL(compressInt(i,bitCount),18u);
    }

    const uint64_t testval(123039843249);
    const uint64_t expect(128849018879);
    BOOST_REQUIRE_EQUAL(compressInt(testval,bitCount),expect);

    // example in function doc:
    BOOST_REQUIRE_EQUAL(compressInt(67u,bitCount),72u);
}

BOOST_AUTO_TEST_CASE( test_bias )
{
    static const unsigned bitCount(2);
    static const unsigned N(1024);
    static const double epsilon(0.000001);

    double sum(0);
    for (unsigned i(0); i<N; ++i)
    {
        sum += compressInt(i,bitCount);
    }
    sum /= N;
    const double expect((N-1)/2.);
    BOOST_REQUIRE_CLOSE(sum,expect,epsilon);
}

BOOST_AUTO_TEST_SUITE_END()

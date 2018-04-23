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

#include "testConfig.h"

#include "htsapi/samtools_fasta_util.hh"
#include "htsapi/vcf_record_util.hh"

#include "common/Exceptions.hh"

#include "boost/test/unit_test.hpp"


BOOST_AUTO_TEST_SUITE( test_vcf_record_util )

static
const char*
getTestVcfPath()
{
    static const std::string testPath(std::string(TEST_DATA_PATH) + "/vcf_record_util_test.vcf.gz");
    return testPath.c_str();
}

static
const char*
getTestFastaPath()
{
    static const std::string testPath(std::string(TEST_DATA_PATH) + "/vcf_record_util_test.fa");
    return testPath.c_str();
}


BOOST_AUTO_TEST_CASE( test_vcf_record_util )
{
    const std::string region("chrP:1-100");

    reference_contig_segment ref;
    get_region_seq(getTestFastaPath(), region, ref.seq());

    vcf_streamer vcfs(getTestVcfPath(),region.c_str());

    const vcf_record* vptr(nullptr);

    // test that 3 VCF records exist and match the reference:
    for (unsigned testIndex(0); testIndex<3; ++testIndex)
    {
        BOOST_REQUIRE(vcfs.next());

        vptr = vcfs.get_record_ptr();
        BOOST_REQUIRE(vptr != nullptr);
        BOOST_REQUIRE(isExpectedVcfReference(ref, *vptr));
        BOOST_REQUIRE_NO_THROW(assertExpectedVcfReference(ref, vcfs));
    }

    // test that 3 VCF records exist and do not match the reference:
    for (unsigned testIndex(0); testIndex<3; ++testIndex)
    {
        BOOST_REQUIRE(vcfs.next());

        vptr = vcfs.get_record_ptr();
        BOOST_REQUIRE(vptr != nullptr);
        BOOST_REQUIRE(! isExpectedVcfReference(ref, *vptr));
        BOOST_REQUIRE_THROW(assertExpectedVcfReference(ref, vcfs), illumina::common::GeneralException);
    }
}

BOOST_AUTO_TEST_SUITE_END()

// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Rumovsky
// Copyright (c) 2013-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

///
/// \author Ole Schulz-Trieglaff
///

#include "GreedyAssembler.hh"

#include "boost/test/unit_test.hpp"

#include <iostream>

BOOST_AUTO_TEST_SUITE( test_GreedyAssembler )

BOOST_AUTO_TEST_CASE( test_runGreedyAssembler )
{
    Assembly as;
    AssemblyReadOutput readInfo;
    AssemblyReadInput reads;
    GreedyAssemblerOptions opt;

    opt.minWordLength         = 5;
    opt.maxWordLength         = 5;
    opt.minContigLength       = 0;
    opt.maxAssemblyIterations = 1;
    opt.minSeedReads          = 1;

    reads.push_back(std::make_pair<int,std::string>(1,"ACACTTTT"));
    reads.push_back(std::make_pair<int,std::string>(2,"TTTTCCAC"));

    runGreedyAssembler(opt,reads,readInfo,as);

    BOOST_REQUIRE_EQUAL(as.size(),1);
    // reverse complement of ACACTTTTTCCAC
    BOOST_REQUIRE_EQUAL(as[0].seq,"GTGGAAAAGTGT");
}

BOOST_AUTO_TEST_SUITE_END()


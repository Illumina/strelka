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

#include "DfsAssembler.hh"

#include "boost/test/unit_test.hpp"

BOOST_AUTO_TEST_SUITE( test_DfsAssembler )


BOOST_AUTO_TEST_CASE( test_runDfsAssembler )
{
    Assembly as;
    AssemblyReadOutput readInfo;
    AssemblyReadInput reads;
    DfsAssemblerOptions opt;

    opt.minWordLength   = 5;
    opt.maxWordLength   = 5;
    opt.minContigLength = 0;
    reads.push_back(std::make_pair<int,std::string>(1,"AAAATTTT"));
    reads.push_back(std::make_pair<int,std::string>(2,"TTTTCCCC"));
    reads.push_back(std::make_pair<int,std::string>(3,"CCCCGGGG"));
    reads.push_back(std::make_pair<int,std::string>(4,"CCCCACAC"));
    reads.push_back(std::make_pair<int,std::string>(5,"CACACGAG"));

    std::string refSeq("AAAATTTTCCCCGGGG");

    runDfsAssembler(opt,reads,refSeq,readInfo,as);

    //for (auto c: as)
    //{
    //    std::cout << "contig=" << c.seq << std::endl;
    //}
    BOOST_REQUIRE_EQUAL(as.size(),2);
    BOOST_REQUIRE_EQUAL(as[0].seq,"AAAATTTTCCCCACACGAG");
    BOOST_REQUIRE_EQUAL(as[1].seq,"AAAATTTTCCCCGGGG");
}

BOOST_AUTO_TEST_CASE( test_runDfsAssembler2 )
{
    Assembly as;
    AssemblyReadOutput readInfo;
    AssemblyReadInput reads;
    DfsAssemblerOptions opt;

    opt.minWordLength   = 5;
    opt.maxWordLength   = 5;
    opt.minContigLength = 0;
    reads.push_back(std::make_pair<int,std::string>(1,"AAAATTTT"));
    reads.push_back(std::make_pair<int,std::string>(2,"TTTTCCCC"));
    reads.push_back(std::make_pair<int,std::string>(3,"CCCCTTTT"));
    reads.push_back(std::make_pair<int,std::string>(4,"TTTTGGGG"));

    std::string refSeq("AAAATTTTCCCCTTTTGGGG");

    runDfsAssembler(opt,reads,refSeq,readInfo,as);

    //for (auto c: as)
    //{
    //    std::cout << "contig=" << c.seq << std::endl;
    //}
    BOOST_REQUIRE_EQUAL(as.size(),1);
    BOOST_REQUIRE_EQUAL(as[0].seq,"AAAATTTTCCCCTTTTGGGG");
}


BOOST_AUTO_TEST_SUITE_END()


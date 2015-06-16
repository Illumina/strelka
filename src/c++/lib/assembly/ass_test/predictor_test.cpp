// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Starka
// Copyright (c) 2009-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

#include "boost/test/unit_test.hpp"
#include "assembly/predictor_bed.hh"

#include <cstdio>

#include "htsapi/hts_streamer.cpp"
#include "htsapi/bed_streamer.cpp"

BOOST_AUTO_TEST_SUITE( assembly_predictor )

BOOST_AUTO_TEST_CASE( test_bed_predictor )
{
    char buffer [L_tmpnam];
    tmpnam (buffer);
    
    FILE * fp = fopen(buffer, "w");
    fprintf(fp, "chr1\t1000\t2000\n");
    fprintf(fp, "chr1\t3000\t4000\n");
    fprintf(fp, "chr1\t5000\t6000\n");
    fprintf(fp, "chr1\t7000\t8000\n");
    fclose(fp);

    predictor_bed bp(buffer, "chr1");

    remove(tmpnam);
}

BOOST_AUTO_TEST_SUITE_END()


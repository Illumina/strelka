// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

#include "boost/test/unit_test.hpp"

#include "cmodel.cpp"


BOOST_AUTO_TEST_SUITE( cmodel )



BOOST_AUTO_TEST_CASE( test_cmodel_qscore ) {

    BOOST_CHECK_EQUAL(prior_adjustment(std::log(2./1.),1.),5);
    BOOST_CHECK_EQUAL(prior_adjustment(std::log(1./2.),1.),2);

    BOOST_CHECK_EQUAL(prior_adjustment(std::log(2000./1.),1.),33);
    BOOST_CHECK_EQUAL(prior_adjustment(std::log(1./2000.),1.),0);

//    BOOST_CHECK_EQUAL(prior_adjustment(std::log(1./0.),1.),40);
}


BOOST_AUTO_TEST_SUITE_END()

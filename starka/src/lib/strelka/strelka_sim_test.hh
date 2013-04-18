// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2012 Illumina, Inc.
//
// This software is covered by the "Illumina Genome Analyzer Software
// License Agreement" and the "Illumina Source Code License Agreement",
// and certain third party copyright/licenses, and any user of this
// source file is bound by the terms therein (see accompanying files
// Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
// Illumina_Source_Code_License_Agreement.pdf and third party
// copyright/license notices).
//
//

/// \file
///
/// \author Chris Saunders
///
/// note coding convention for all ranges '_pos fields' is:
/// XXX_begin_pos is zero-indexed position at the begining of the range
/// XXX_end_pos is zero-index position 1 step after the end of the range
///

#ifndef STRELKA_SIM_TEST_
#define STRELKA_SIM_TEST_


enum sim_mode_t {
    SIM_RANDOM,
    SIM_NOISE,
    SIM_SOMATIC,
    SIM_GERMLINE,
    SIM_REF
};


// the internal strelka simulator is not designed to be a fully
// fleshed out tool... this is just the minimum we need to do
// fast site simulations...
//
struct strelka_site_sim_options {

    strelka_site_sim_options()
        : total_sites(1000000)
        , ncov(35)
        , tcov(35)
        , qval_file("default_data/qval_distro")
        , ssnv_prior(0.000001)
        , tumor_purity(1.)
        , seed(1)
        , mode(SIM_RANDOM)
    {}

    unsigned total_sites;
    unsigned ncov;
    unsigned tcov;
    std::string qval_file;
    std::string oracle_file;

    double ssnv_prior;
    double tumor_purity;
    uint32_t seed;
    sim_mode_t mode;
};


// completely self-contained version of pile simulator:
//
void
strelka_site_sim(strelka_options& opt,
                 strelka_site_sim_options& sim_opt);

#endif

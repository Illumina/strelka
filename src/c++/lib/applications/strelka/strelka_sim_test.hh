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
/// \author Chris Saunders
///

#pragma once


enum sim_mode_t
{
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
struct strelka_site_sim_options
{

    strelka_site_sim_options()
        : total_sites(1000000)
        , ncov(35)
        , tcov(35)
        , ssnv_prior(0.000001)
        , normal_purity(1.)
        , tumor_purity(1.)
        , seed(1)
        , mode(SIM_RANDOM)
        , is_somatic_gvcf(false)
    {}

    unsigned total_sites;
    unsigned ncov;
    unsigned tcov;
    std::string qval_file;
    std::string oracle_file;

    double ssnv_prior;
    double normal_purity;
    double tumor_purity;
    uint32_t seed;
    sim_mode_t mode;

    bool is_somatic_gvcf;
};


// completely self-contained version of pile simulator:
//
void
strelka_site_sim(strelka_options& opt,
                 strelka_site_sim_options& sim_opt);

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

///
/// \author Chris Saunders
///

#pragma once

#include "starling_shared.hh"

#include <string>


enum sim_mode_t
{
    SIM_RANDOM,
    SIM_NOISE,
    SIM_GERMLINE,
    SIM_REF
};


// the internal starling simulator is not designed to be a fully
// fleshed out tool... this is just the minimum we need to do
// fast site simulations...
//
struct starling_site_sim_options
{
    unsigned total_sites = 1000000;
    unsigned coverage = 35;
    std::string qval_file;
    std::string oracle_file;

    uint32_t seed = 1;
    sim_mode_t mode = SIM_RANDOM;
    bool is_het_only = false;
    bool is_exact_cov = false;
};


void
starling_site_sim(
    starling_options& opt,
    starling_site_sim_options& sim_opt);

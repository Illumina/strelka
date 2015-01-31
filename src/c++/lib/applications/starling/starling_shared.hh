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

///
/// \author Chris Saunders
///

#pragma once

//#include "gvcf_options.hh"
#include "starling_common/starling_base_shared.hh"



struct starling_options : public starling_base_options
{
    typedef starling_options base_t;

    starling_options() {}

    bool
    is_bsnp_diploid() const override { return true; }

    //gvcf_options gvcf;
};


// data deterministically derived from the input options:
//
struct starling_deriv_options : public starling_base_deriv_options
{
    typedef starling_base_deriv_options base_t;

    starling_deriv_options(
        const starling_options& opt,
        const reference_contig_segment& ref)
        : base_t(opt,ref)
    {}

//    gvcf_deriv_options gvcf;
};


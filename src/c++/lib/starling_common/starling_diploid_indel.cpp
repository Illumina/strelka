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

#include "starling_diploid_indel.hh"

#include <iostream>


void
starling_diploid_indel_core::
dump(std::ostream& os) const
{
   os << "ploidy=" << ploidy
        << ",max_gt=" << max_gt
        << ",indel_qphred=" << indel_qphred
        << ",max_gt_qphred=" << max_gt_qphred
        << ",max_gt_poly=" << max_gt_poly
        << ",max_gt_poly_qphred=" << max_gt_poly_qphred
        << ",is_forced_output=" << is_forced_output
        << ",is_zero_coverage=" << is_zero_coverage;
}



void
starling_diploid_indel::
dump(std::ostream& os) const
{
    starling_diploid_indel_core::dump(os);
    os << ",pprob=[";

    bool isFirst(true);
    for (const auto p : pprob)
    {
        if (! isFirst) os << ',';
        os << p;
        isFirst=false;
    }
    os << ']';
}



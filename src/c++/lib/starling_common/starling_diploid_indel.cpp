// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2016 Illumina, Inc.
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

#include "starling_diploid_indel.hh"

#include <iostream>



const int starling_diploid_indel_core::maxQ = 999;



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



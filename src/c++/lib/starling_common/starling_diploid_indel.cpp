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



const int GermlineDiploidIndelSimpleGenotypeInfoCore::maxQ = 999;



void
GermlineDiploidIndelSimpleGenotypeInfoCore::
dump(std::ostream& os) const
{
    os << "max_gt=" << max_gt
       << ",max_gt_qphred=" << max_gt_qphred
       << ",max_gt_poly=" << max_gt_poly
       << ",max_gt_poly_qphred=" << max_gt_poly_qphred;
}



void
starling_diploid_indel::
dump(std::ostream& os) const
{
    GermlineDiploidIndelSimpleGenotypeInfoCore::dump(os);
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



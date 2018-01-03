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

#include "blt_common/snp_pos_info.hh"
#include "blt_util/monogt.hh"

#include <boost/utility.hpp>

#include <iosfwd>


struct monoploid_genotype : private boost::noncopyable
{
    monoploid_genotype() : is_snp(false), ref_gt(0), max_gt(0), max2_gt(0)
    {
        static const double p(1./static_cast<double>(MONOGT::SIZE));
        for (unsigned i(0); i<MONOGT::SIZE; ++i) pprob[i] = p;
    }

    bool is_snp;
    unsigned ref_gt;
    unsigned max_gt;
    unsigned max2_gt;
    double pprob[MONOGT::SIZE];
};

std::ostream& operator<<(std::ostream& os,monoploid_genotype& dp);



/// \brief call a snp @ pos by calculating the posterior probability
/// of all possible genotypes for a monoploid individual.
///
void
position_snp_call_pprob_monogt(const double theta,
                               const snp_pos_info& pi,
                               monoploid_genotype& dgt);

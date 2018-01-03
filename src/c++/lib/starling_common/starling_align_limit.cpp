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
///
/// \author Chris Saunders
///

#include "starling_common/starling_align_limit.hh"

#include "boost/math/special_functions/binomial.hpp"

#include <cmath>



#if 0
// Calculate the maximum number of candidate alignments for a given
// set of togglable indels and a toggle depth.
//
static
float
max_candidate_alignment_count(const unsigned n_indel,
                              const unsigned n_toggle)
{
    const float n(n_indel);

    float sum(1);
    for (unsigned i(0); i<n_toggle; ++i)
    {
        const float k(i+1);
        sum += std::pow(2.,k)*boost::math::binomial_coefficient<float>(n,k);
    }
    return sum;
}
#endif



// Calculate the number of indel toggles which cannot exceed the
// maximum candidate alignment count.
//
static
unsigned
max_candidate_alignment_toggle(const unsigned n_indel,
                               const unsigned max_alignments)
{
    const float max(max_alignments);

    float sum(1.);
    for (unsigned i(0); i<n_indel; ++i)
    {
        const unsigned k(i+1);
        sum += std::pow(static_cast<float>(2),static_cast<float>(k))*boost::math::binomial_coefficient<float>(n_indel,k);
        if (sum>max) return i;
    }
    return n_indel;
}



starling_align_limit::
starling_align_limit(const unsigned mac)
{
    static const unsigned max_indels(100);
    for (unsigned i(0); i<max_indels; ++i)
    {
        const unsigned mt(max_candidate_alignment_toggle(i,mac));
        if ((i>1) && (mt<2)) break;
        _max_toggle.push_back(mt);
    }
}

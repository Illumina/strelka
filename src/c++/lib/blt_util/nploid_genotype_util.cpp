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
#include "blt_util/nploid_genotype_util.hh"

#include <boost/math/special_functions/binomial.hpp>



unsigned
nploid_gtype_size(const unsigned ploidy)
{

    static const unsigned NB1(N_BASE-1);
    const double r(boost::math::binomial_coefficient<double>(NB1+ploidy,NB1));
    assert(r>=0);
    return static_cast<unsigned>(r);
}



nploid_info::
nploid_info(const unsigned init_ploidy)
    : _ploidy(init_ploidy), _echunk(1./static_cast<double>(init_ploidy))
{

    gtype_info gi;
    gi.label.resize(_ploidy);
    ginfo_init(_ploidy,0,0,gi);

    unsigned init_index(0);
    unsigned* init_ref(_ref_gtype);
    ref_gtype_init(_ploidy,0,0,true,init_ref,init_index);
}



// recursively derive labels and base frequencies for each genotype
void
nploid_info::
ginfo_init(const unsigned init_ploidy,
           const unsigned pli,
           const unsigned init_i,
           gtype_info& gi)
{

    for (unsigned i(init_i); i<N_BASE; ++i)
    {
        gi.efreq_levels[i] += 1;
        gi.label[pli]=id_to_base(i);
        if ((pli+1)<init_ploidy)
        {
            ginfo_init(init_ploidy,pli+1,i,gi);
        }
        else
        {
            _ginfo.push_back(gi);
        }
        gi.efreq_levels[i] -= 1;
    }
}



// build an index of which genotype corresponds to the reference state,
// for each possible reference base:
void
nploid_info::
ref_gtype_init(const unsigned init_ploidy,
               const unsigned pli,
               const unsigned init_i,
               const bool is_hom,
               unsigned*& ref,
               unsigned& index)
{

    for (unsigned i(init_i); i<N_BASE; ++i)
    {
        const bool is_iter_hom(pli==0 || (is_hom && (i==init_i)));
        if ((pli+1)<init_ploidy)
        {
            ref_gtype_init(init_ploidy,pli+1,i,is_iter_hom,ref,index);
        }
        else
        {
            if (is_iter_hom) *ref++ = index;
            ++index;
        }
    }
}


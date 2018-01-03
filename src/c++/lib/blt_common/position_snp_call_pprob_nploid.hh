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

#include "blt_common/snp_pos_info.hh"

#include "blt_util/nploid_genotype_util.hh"

#include <iosfwd>


struct nploid_genotype
{

    nploid_genotype(const nploid_info& ni)
        : is_snp(false), ref_gt(0), max_gt(0), max2_gt(0),
          pprob(ni.gtype_size()), ploidy(ni.ploidy())
    {
        const unsigned n_gt(ni.gtype_size());
        static const double p(1./static_cast<double>(n_gt));
        for (unsigned i(0); i<n_gt; ++i) pprob[i] = p;
    }

    bool is_snp;
    unsigned ref_gt;
    unsigned max_gt;
    unsigned max2_gt;
    std::vector<double> pprob;
    unsigned ploidy;
};


void
nploid_write(const nploid_info& ninfo,
             const nploid_genotype& ng,
             std::ostream& os);


/// \brief call a snp @ pos by calculating the posterior probability
/// of all possible genotypes for an n-ploid individual sample, ploidy
/// is specified by nploid_info.
///
void
position_snp_call_pprob_nploid(const double snp_prob,
                               const snp_pos_info& pi,
                               const nploid_info& ninfo,
                               nploid_genotype& ngt);

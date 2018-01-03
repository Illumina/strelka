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

#include "blt_common/blt_shared.hh"

#include <boost/utility.hpp>

#include <iosfwd>


struct lrt_snp_call : private boost::noncopyable
{

    lrt_snp_call()
        : is_snp(false),
          null_loghood(0),
          min_test_loghood(0),
          min_loghood(0),
          snp_prob(0)
    {
        for (unsigned i(0); i<N_BASE; ++i) allele_freq[i]=0.;
    }

    bool is_snp;
    double null_loghood;
    double min_test_loghood; // lhood from the single-parameter test used to determine significance
    double min_loghood;      // lhood from the three-parameter minimization used to find all ML frequencies
    double snp_prob;
    double allele_freq[N_BASE]; // ML frequencies taking Qvalues into account.
};

std::ostream& operator<<(std::ostream& os,
                         const lrt_snp_call& sc);



/// \brief call a snp @ pos by likelihood ratio test
///
/// just returns yes/no for a snp right now, this will have to return
/// some more useful info eventually...
///
void
position_snp_call_lrt(const double alpha,
                      const snp_pos_info& pi,
                      lrt_snp_call& sc);

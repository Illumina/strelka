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

#include "blt_common/blt_shared.hh"
#include "blt_common/nonref_test_call.hh"
#include "blt_common/snp_pos_info.hh"

#include <iosfwd>


namespace NRTEST
{
enum index_t
{
    REF,
    NONREF,
    SIZE
};

inline
const char*
label(const index_t i)
{
    switch (i)
    {
    case REF:
        return "ref";
    case NONREF:
        return "nonref";
    default:
        return "xxx";
    }
}
}


/// \brief Call a snp at a position under the assumption that any
/// combination of non-reference requencies could occur. When a snp is
/// found, optionally also provide the allele MLEs.
///
void
position_nonref_test(const snp_pos_info& pi,
                     const double nonref_site_prob,
                     const double min_nonref_freq,
                     const bool is_mle_freq,
                     nonref_test_call& nrc);


void
write_nonref_test(const blt_options& opt,
                  const snp_pos_info& pi,
                  const nonref_test_call& nrc,
                  std::ostream& os);

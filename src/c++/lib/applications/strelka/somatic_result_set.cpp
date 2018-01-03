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
/// \author Chris Saunders, Sangtae Kim
///


#include "somatic_result_set.hh"

#include "strelka_digt_states.hh"
#include "somatic_call_shared.hh"

#include <iostream>



std::ostream&
operator<<(std::ostream& os, const result_set& rs)
{
    os << "rs: ntype: " << NTYPE::label(rs.ntype) << " qphred: " << rs.qphred << " n_qphred " << rs.from_ntype_qphred << "\n";
    os << " max_gt: ";
    DDIGT::write_indel_state(static_cast<DDIGT::index_t>(rs.max_gt),os);
    return os;
}

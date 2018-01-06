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

#include "IndelKey.hh"
#include "common/Exceptions.hh"

#include <iostream>
#include <sstream>



void
IndelKey::
validate() const
{
    using namespace illumina::common;

    // insertSequence/deletionlength is not intended for breakpoint types:
    if (is_breakpoint())
    {
        if ((delete_length() != 0) or (insert_length() != 0))
        {
            std::ostringstream oss;
            oss << "Invalid allele type -- breakpoint also has insertion/deletion defined '" << (*this) << "'";
            BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
        }
    }

#if 0
    // an indel with zero insert and delete length is assumed to be an error:
    if (type == INDEL::INDEL)
    {
        if ((delete_length() == 0) and (insert_length() == 0))
        {
            std::ostringstream oss;
            oss << "ERROR: invalid allele type -- indel with no insertion/deletion variant defined '" << (*this) << "'\n";
            BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
        }
    }
#endif
}

//bool
//IndelKey::
//isOrthogonal(IndelKey other) const
//{
//    // indels belonging to different active regions are not orthogonal
//    if (activeRegionId == 0 || activeRegionId != other.activeRegionId) return false;
//
//
//    return false;
//}

std::ostream&
operator<<(
    std::ostream& os,
    const IndelKey& indelKey)
{
    os << "INDEL pos: " << indelKey.pos
       << " type: " << INDEL::get_index_label(indelKey.type)
       << " deleteLength: " << indelKey.delete_length()
       << " insertSeq: " << indelKey.insert_seq() << "\n";
    return os;
}



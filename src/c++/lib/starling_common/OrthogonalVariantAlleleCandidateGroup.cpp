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


#include "OrthogonalVariantAlleleCandidateGroup.hh"

#include <iostream>



known_pos_range
OrthogonalVariantAlleleCandidateGroup::
getReferenceRange() const
{
    known_pos_range pr(0,0);
    const unsigned alleleCount(size());
    assert(alleleCount!=0);

    for (unsigned alleleIndex(0); alleleIndex<alleleCount; alleleIndex++)
    {
        const IndelKey& indelKey(key(alleleIndex));
        const known_pos_range indelRange(indelKey.pos, indelKey.right_pos());
        if (alleleIndex==0)
        {
            pr = indelRange;
        }
        else
        {
            pr.merge_range(indelRange);
        }
    }
    return pr;
}



std::ostream&
operator<<(
    std::ostream& os, const OrthogonalVariantAlleleCandidateGroup& group)
{
    const unsigned altAlleleCount(group.size());
    os << "orthogonalAlleleGroup size: " << altAlleleCount << "  keys: \n";
    for (unsigned altAlleleIndex(0); altAlleleIndex < altAlleleCount; ++altAlleleIndex)
    {
        os << group.key(altAlleleIndex) << "\n";
    }
    return os;
}

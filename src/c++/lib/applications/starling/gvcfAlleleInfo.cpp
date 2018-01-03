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


#include "gvcfAlleleInfo.hh"
#include "blt_util/math_util.hh"
#include "common/Exceptions.hh"

#include "boost/math/distributions/binomial.hpp"
#include "boost/math/distributions.hpp"
#include "rnaVariantEmpiricalScoringFeatures.hh"

#include <iostream>
#include <map>
#include <sstream>
#include <typeinfo>


/// remove numPrefixToRemove matches from a cigar match segment of length matchSegmentLength
static
void
removeCommonPrefixFromCigar(unsigned& matchSegmentLength, unsigned& numPrefixToRemove, bool isFirstSegment = false)
{
    if (isFirstSegment)
        --matchSegmentLength;   // save the initial match

    if (matchSegmentLength >= numPrefixToRemove)
    {
        matchSegmentLength -= numPrefixToRemove;
        numPrefixToRemove = 0;
    }
    else
    {
        numPrefixToRemove -= matchSegmentLength;
        matchSegmentLength = 0;
    }

    if (isFirstSegment)
        ++matchSegmentLength;   // restore the initial match
}



void
setIndelAlleleCigar(
    unsigned lead,
    unsigned trail,
    const unsigned commonPrefixLength,
    const IndelKey& indelKey,
    ALIGNPATH::path_t& cigar)
{
    using namespace ALIGNPATH;

    unsigned numPrefixToRemove(commonPrefixLength);
    removeCommonPrefixFromCigar(lead, numPrefixToRemove, true);

    cigar.clear();
    if (lead)
    {
        cigar.push_back(path_segment(MATCH,lead));
    }
    if (indelKey.delete_length())
    {
        cigar.push_back(path_segment(DELETE, indelKey.delete_length()));
    }
    unsigned insertLength(indelKey.insert_length());
    if (numPrefixToRemove)
        removeCommonPrefixFromCigar(insertLength, numPrefixToRemove);

    if (insertLength)
        cigar.push_back(path_segment(INSERT, insertLength));

    if (numPrefixToRemove)
        removeCommonPrefixFromCigar(trail, numPrefixToRemove);

    // numPrefixToRemove > 0 means that commonPrefixLength is incorrect
    assert(numPrefixToRemove == 0);

    if (trail)
    {
        cigar.push_back(path_segment(MATCH,trail));
    }
}



std::ostream&
operator<<(
    std::ostream& os,
    const GermlineVariantAlleleInfo& allele)
{
    os << "isForcedOutput: " << allele.isForcedOutput;
    return os;
}



std::ostream&
operator<<(
    std::ostream& os,
    const GermlineSiteAlleleInfo& allele)
{
    os << static_cast<GermlineVariantAlleleInfo>(allele) << '\n';
    os << "base: " << id_to_base(allele.baseIndex);
    return os;
}



std::ostream&
operator<<(
    std::ostream& os,
    const GermlineIndelAlleleInfo& allele)
{
    os << static_cast<GermlineVariantAlleleInfo>(allele) << '\n';
    os << "IndelKey: " << allele.indelKey << "\n";
    os << "AlleleReportInfo: " << allele.indelReportInfo << "\n";
    return os;
}

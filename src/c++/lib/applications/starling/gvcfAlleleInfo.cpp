// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2016 Illumina, Inc.
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





void
setIndelAlleleCigar(
    const unsigned lead,
    const unsigned trail,
    const IndelKey& indelKey,
    ALIGNPATH::path_t& cigar)
{
    using namespace ALIGNPATH;

    cigar.clear();
    if (lead)
    {
        cigar.push_back(path_segment(MATCH,lead));
    }
    if (indelKey.delete_length())
    {
        cigar.push_back(path_segment(DELETE, indelKey.delete_length()));
    }
    if (indelKey.insert_length())
    {
        cigar.push_back(path_segment(INSERT, indelKey.insert_length()));
    }
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
    os << "AleleReportInfo: " << allele.indelReportInfo << "\n";
    return os;
}

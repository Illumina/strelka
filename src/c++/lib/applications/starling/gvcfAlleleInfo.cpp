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
GermlineIndelAlleleInfo::
set_hap_cigar(
    const unsigned lead,
    const unsigned trail)
{
    using namespace ALIGNPATH;

    cigar.clear();
    if (lead)
    {
        cigar.push_back(path_segment(MATCH,lead));
    }
    if (_indelKey.delete_length())
    {
        cigar.push_back(path_segment(DELETE,_indelKey.delete_length()));
    }
    if (_indelKey.insert_length())
    {
        cigar.push_back(path_segment(INSERT,_indelKey.insert_length()));
    }
    if (trail)
    {
        cigar.push_back(path_segment(MATCH,trail));
    }
}



std::ostream&
operator<<(std::ostream& os,
           const GermlineVariantAlleleInfo& shmod)
{
    os << "gqx: " << shmod.gqx
       << " gq: " << shmod.gq;
    if (typeid(shmod) == typeid(GermlineDiploidIndelAlleleInfo))
    {
        auto imod = dynamic_cast<const GermlineDiploidIndelAlleleInfo&>(shmod);

        os << " max_gt: " << DIGT::label(imod.max_gt);
    }

    return os;
}

std::ostream&
operator<<(std::ostream& os,
           const GermlineDiploidSiteAlleleInfo& smod)
{
    os << static_cast<GermlineVariantAlleleInfo>(smod) << '\n';

    os << "is_unknown: " << smod.is_unknown;
    os << " is_covered: " << smod.is_covered;
    os << " is_used_coverage: " << smod.is_used_covered;
    os << " is_zero_ploidy: " << smod.is_zero_ploidy;

    if (smod.modified_gt != MODIFIED_SITE_GT::NONE)
    {
        os << " modgt: " << MODIFIED_SITE_GT::get_label(smod.modified_gt);
    }

    return os;
}



std::ostream&
operator<<(
    std::ostream& os,
    const GermlineIndelAlleleInfo& shi)
{
    os << static_cast<GermlineVariantAlleleInfo>(shi) << '\n';

    os << "IndelKey: " << shi._indelKey << "\n";
    //os << "indel_data: " << shi._id << "\n";
    os << "indel_report_info: " << shi._indelReportInfo << "\n";
    os << "indel_sample_info: " << shi._indelSampleReportInfo << "\n";
    os << "cigar: " << shi.cigar << "\n";

    return os;
}



std::ostream&
operator<<(
    std::ostream& os,
    const GermlineDiploidIndelAlleleInfo& dic)
{
    os << static_cast<GermlineIndelAlleleInfo>(dic) << '\n';

    dic._dindel.dump(os);

    os << " max_gt: " << dic.max_gt << "\n";

    return os;
}

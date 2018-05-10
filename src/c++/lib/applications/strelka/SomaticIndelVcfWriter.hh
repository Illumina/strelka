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

#include "somatic_result_set.hh"
#include "strelka_shared.hh"

#include "starling_common/AlleleReportInfo.hh"
#include "../../starling_common/LocalRegionStats.hh"

#include <array>
#include <map>


/// store vcf indel record info to write out later:
struct SomaticIndelVcfInfo
{
    somatic_indel_call sindel;
    AlleleReportInfo indelReportInfo;
    std::string vcf_ref_seq; ///< vcf REF string
    std::string vcf_indel_seq; ///< vcf ALT string
    std::array<AlleleSampleReportInfo,2> nisri;
    std::array<AlleleSampleReportInfo,2> tisri;
};


/// delay writing indel vcf record until window data is computed:
struct SomaticIndelVcfWriter
{
    SomaticIndelVcfWriter(
        const strelka_options& opt,
        const strelka_deriv_options& dopt,
        std::ostream* osptr) :
        _opt(opt),
        _dopt(dopt),
        _osptr(osptr) {}

    /// return true if indel information is cached for this position
    bool
    testPos(
        const pos_t pos) const
    {
        return (_data.count(pos) != 0);
    }

    /// return true if no indel information is cached
    bool
    empty() const
    {
        return (_data.empty());
    }

    void
    clear()
    {
        _data.clear();
    }

    /// store an indel call
    void
    cacheIndel(
        const pos_t pos,
        const SomaticIndelVcfInfo& siInfo);

    /// add final information required
    ///
    /// \param[in] maxChromDepth max expected normal sample depth for this chromosome
    void
    addIndelWindowData(
        const std::string& chromName,
        const pos_t pos,
        const LocalRegionStats& wasNormal,
        const LocalRegionStats& wasTumor,
        const double maxChromDepth);

private:
    const strelka_options& _opt;
    const strelka_deriv_options& _dopt;
    std::ostream* _osptr;
    std::map<pos_t,std::vector<SomaticIndelVcfInfo>> _data;
};

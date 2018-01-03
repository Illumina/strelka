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

#include "snoise_shared.hh"
#include "snoise_streams.hh"
#include "starling_common/starling_pos_processor_base.hh"


///
///
struct snoise_pos_processor : public starling_pos_processor_base
{
    typedef starling_pos_processor_base base_t;

    snoise_pos_processor(
        const snoise_options& opt,
        const starling_base_deriv_options& dopt,
        const reference_contig_segment& ref,
        const snoise_streams& fileStreams,
        RunStatsManager& statsManager);

    void
    resetRegion(
        const std::string& chromName,
        const known_pos_range2& reportRange)
    {
        base_t::resetRegionBase(chromName, reportRange);
    }

private:
    void
    process_pos_variants_impl(
        const pos_t pos,
        const bool isPosPrecedingReportableRange) override
    {
        if (isPosPrecedingReportableRange) return;
        process_pos_snp_snoise(pos);
    }

    void
    process_pos_snp_snoise(const pos_t pos);

    const snoise_streams& _fileStreams;
};

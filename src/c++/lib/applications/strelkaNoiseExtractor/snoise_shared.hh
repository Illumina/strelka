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

#pragma once

#include "starling_common/starling_base_shared.hh"

struct snoise_options : public starling_base_options
{
    snoise_options()
    {
        // parameter defaults are similar to the germline caller
        // TODO: should this be updated with all of the latest germline caller settings for haplotyping, etc?
        bsnp_ssd_no_mismatch = 0.35;
        bsnp_ssd_one_mismatch = 0.6;
        mismatchDensityFilterMaxMismatchCount = 2;
        mismatchDensityFilterFlankSize = 20;
        is_min_vexp = true;
        min_vexp = 0.25;
    }

    const AlignmentFileOptions&
    getAlignmentFileOptions() const override
    {
        return alignFileOpt;
    }

    AlignmentFileOptions alignFileOpt;

    bool is_skip_header = false;
};

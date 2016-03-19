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

#pragma once

#include "blt_util/chrom_depth_map.hh"
#include "starling_common/starling_base_shared.hh"


struct SequenceErrorCountsOptions : public starling_base_options
{
    typedef SequenceErrorCountsOptions base_t;

    SequenceErrorCountsOptions()
    {
        // manipulate indel candidate defaults, so that EVERYTHING is a candidate:
        is_candidate_indel_signal_test = false;

        // filter candidates based on depth more stringently than we do in germline
        // calling, this helps to partially offset the increase in candidates from
        // turning the signal test off (above)
        max_candidate_indel_depth_factor = 1;

        // set command-line defaults for starling only:
        min_paired_align_score = 20;
    }

    bool
    is_depth_filter() const
    {
        return (! chrom_depth_file.empty());
    }

    std::string chrom_depth_file;
    double max_depth_factor = 3.;

    std::string countsFilename;
};


// data deterministically derived from the input options:
//
struct SequenceErrorCountsDerivOptions : public starling_base_deriv_options
{
    typedef starling_base_deriv_options base_t;

    SequenceErrorCountsDerivOptions(
        const SequenceErrorCountsOptions& opt,
        const reference_contig_segment& ref);

    bool
    is_max_depth() const
    {
        return (! chrom_depth.empty());
    }

    std::string block_label;
    double max_depth = 0;
    double norm_depth = 0;
    cdmap_t chrom_depth;
};


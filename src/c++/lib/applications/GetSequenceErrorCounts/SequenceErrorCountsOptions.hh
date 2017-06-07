//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2017 Illumina, Inc.
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

        // set command-line defaults for error counting only:
        min_mapping_quality = 60;
        min_qscore = 17;
        isBasecallQualAdjustedForMapq = false;

        is_max_win_mismatch = true;
        max_win_mismatch = 1;
        max_win_mismatch_flank_size = 100;

        tier2.is_tier2_min_mapping_quality = true;
        tier2.tier2_min_mapping_quality = 0;

        minDistanceFromReadEdge = 10;
    }

    bool
    is_depth_filter() const
    {
        return (! chrom_depth_file.empty());
    }

    bool is_write_observations() const
    {
        return (! observationsBedFilename.empty());
    }

    const AlignmentFileOptions&
    getAlignmentFileOptions() const override
    {
        return alignFileOpt;
    }

    AlignmentFileOptions alignFileOpt;

    std::string chrom_depth_file;
    double max_depth_factor = 3.;

    //========= input files:
    /// File indicating known variant genotypes
    std::string knownVariantsFile;

    /// File indicating which regions to exclude from counting
    std::vector<std::string> excludedRegionsFileList;

    //======== output files:
    /// Binary output file containing error counts data (required)
    std::string countsFilename;

    /// optional debug output specifying indel location
    std::string observationsBedFilename;

    /// optional evidence count indicating the number of non-empty sites considered during error counting
    std::string nonEmptySiteCountFilename;
};


/// data deterministically derived from the input options:
///
struct SequenceErrorCountsDerivOptions : public starling_base_deriv_options
{
    typedef starling_base_deriv_options base_t;

    explicit
    SequenceErrorCountsDerivOptions(const SequenceErrorCountsOptions& opt);

    bool
    is_max_depth() const
    {
        return (! chrom_depth.empty());
    }

    std::string block_label;
    cdmap_t chrom_depth;
};

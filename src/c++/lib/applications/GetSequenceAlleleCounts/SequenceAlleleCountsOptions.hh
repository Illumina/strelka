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

#include "blt_util/chrom_depth_map.hh"
#include "starling_common/starling_base_shared.hh"


struct SequenceAlleleCountsOptions : public starling_base_options
{
    typedef SequenceAlleleCountsOptions base_t;

    SequenceAlleleCountsOptions()
    {
        // manipulate indel candidate defaults, so that EVERYTHING is a candidate:
        is_candidate_indel_signal_test = false;

        // filter candidates based on depth more stringently than we do in germline
        // calling, this helps to partially offset the increase in candidates from
        // turning the signal test off (above)
        max_candidate_indel_depth_factor = 1;

        // set command-line defaults for error counting only:
        minMappingErrorPhredProb = 60;
        minBasecallErrorPhredProb = 17;
        isBasecallQualAdjustedForMapq = false;

        // custom mmdf for sequence error counting
        mismatchDensityFilterMaxMismatchCount = 1;
        mismatchDensityFilterFlankSize = 100;

        // use tier2 but disable everything besides the lower MAPQ value
        useTier2Evidence = true;
        tier2.mismatchDensityFilterMaxMismatchCount = 2;
        tier2.includeAnomalousReads = false;
        tier2.includeSingletonReads = false;
        tier2.isRandomBaseMatchProb = false;

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
struct SequenceAlleleCountsDerivOptions : public starling_base_deriv_options
{
    typedef starling_base_deriv_options base_t;

    explicit
    SequenceAlleleCountsDerivOptions(const SequenceAlleleCountsOptions& opt);

    bool
    is_max_depth() const
    {
        return (! chrom_depth.empty());
    }

    std::string block_label;
    cdmap_t chrom_depth;
};

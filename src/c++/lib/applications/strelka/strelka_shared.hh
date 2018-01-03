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

#include "blt_util/chrom_depth_map.hh"
#include "calibration/VariantScoringModelServer.hh"
#include "options/TumorNormalAlignmentFileOptions.hh"
#include "starling_common/starling_base_shared.hh"


/// variant call filtration options used only for somatic snvs and indels
///
/// these are isolated from starling gvcf options to avoid conflicts
///
struct somatic_filter_options
{
    bool
    is_depth_filter() const
    {
        return (! chrom_depth_file.empty());
    }

    std::string chrom_depth_file;
    bool is_skip_header = false;
    double max_depth_factor = 3.;
    double snv_max_filtered_basecall_frac = 0.4;
    double snv_max_spanning_deletion_frac = 0.75;
    int snv_min_qss_ref = 15;

    double indelMaxWindowFilteredBasecallFrac = 0.3;
    int sindelQuality_LowerBound = 30;

    unsigned indelRegionFlankSize = 50;
    double snvMinEVS = 20.0;

    /// LowDepth filter is added if the tumor tier 1 depth is below this value
    unsigned minPassedCallDepth = 2u;
};



struct strelka_options : public starling_base_options
{
    typedef starling_base_options base_t;

    strelka_options()
    {
        // custom mmdf defaults for somatic calling
        mismatchDensityFilterMaxMismatchCount = 3;
        mismatchDensityFilterFlankSize = 20;

        // custom noisy read handling for somatic calling
        randomBaseMatchProb = 0.5;

        minBasecallErrorPhredProb = 0;

        useTier2Evidence = true;

        // turn on empirical scoring for strelka only:
        is_compute_somatic_scoring_metrics = true;

        // Set this true to use haplotyping in somatic calling
        isHaplotypingEnabled = false;

        isSomaticCallingMode = true;
    }

    bool is_somatic_snv() const
    {
        return (! somatic_snv_filename.empty());
    }

    bool is_somatic_indel() const
    {
        return (! somatic_indel_filename.empty());
    }

    bool
    is_somatic_callable() const
    {
        return (! somatic_callable_filename.empty());
    }

    bool
    isUseSomaticSNVScoring() const
    {
        return (! somatic_snv_scoring_model_filename.empty());
    }

    bool
    isUseSomaticIndelScoring() const
    {
        return (! somatic_indel_scoring_model_filename.empty());
    }

    const AlignmentFileOptions&
    getAlignmentFileOptions() const override
    {
        return alignFileOpt;
    }

    TumorNormalAlignmentFileOptions alignFileOpt;

    /// Expected rate of somatic SNVs
    double somatic_snv_rate = 0.000001;
    std::string somatic_snv_filename;

    /// Expected rate of somatic indels
    double somatic_indel_rate = 0.000001;
    std::string somatic_indel_filename;

    double shared_site_error_rate = 0.000005;
    double shared_site_error_strand_bias_fraction = 0.5;

    /// \brief Controls the expectation that an indel observed in the tumor may be a noise allele shared with the normal
    /// sample.
    ///
    /// A higher value means a lower expected shared error rate. Note this is the indel analog of shared_site_error_rate
    /// used for sites. Indels use a different scheme than sites to effectively allow the shared error rate to vary as
    /// a function of each indel's error rate.
    double shared_indel_error_factor = 1.65;

    /// \brief Highest expected contamination level of tumor cells into the normal sample.
    ///
    /// Such contamination is typical for normal controls from liquid and late-stage cancers. The model integrates over
    /// a contamination range from 0 up to this value. The values can be adjusted separately for SNVs and indels,
    /// although they are expected to typically be set to the same value.
    ///@{
    double ssnv_contam_tolerance = 0.15;
    double indel_contam_tolerance = 0.15;
    ///@}

    /// \brief The minimum indel breakpoint flanking distance for a read to be used as indel evidence
    /// in the normal sample.
    ///
    /// A lower flank requirement is used for normal sample reads during somatic variant calling to ensure that all
    /// evidence potentially used against a somatic call in the normal is available.
    int normal_sample_min_read_bp_flank = 1;

    std::string somatic_callable_filename;

    /// \brief Variants in these vcfs are used to indicate known systematic low-freqeuncy noise.
    std::vector<std::string> noise_vcf;

    somatic_filter_options sfilter;

    /// somatic scoring models:
    std::string somatic_snv_scoring_model_filename;
    std::string somatic_indel_scoring_model_filename;

    const std::string SomaticEVSVcfInfoTag = "SomaticEVS";
};



/// somatic filter options computed after user input is finished:
///
struct somatic_filter_deriv_options
{
    bool
    is_max_depth() const
    {
        return (! chrom_depth.empty());
    }

    cdmap_t chrom_depth;

    /// one somatic indel feature requires a delayed computation (it is based on a window around the variant call,
    /// this is accomplished by using a "post-call" stage in the pos-processor stage-manager scheme, here we
    /// store the stage number associated with this post-call filtration step:
    unsigned indelRegionStage = 0;
};



struct somatic_snv_caller_strand_grid;
struct somatic_indel_caller_grid;


/// data deterministically derived from the input options or read in from model files, etc.
///
struct strelka_deriv_options : public starling_base_deriv_options
{
    typedef starling_base_deriv_options base_t;

    explicit
    strelka_deriv_options(const strelka_options& opt);

    ~strelka_deriv_options();

    const somatic_snv_caller_strand_grid&
    sscaller_strand_grid() const
    {
        return *(_sscaller_strand_grid.get());
    }

    const somatic_indel_caller_grid&
    sicaller_grid() const
    {
        return *(_sicaller_grid.get());
    }

/// data:
    somatic_filter_deriv_options sfilter;

    std::unique_ptr<VariantScoringModelServer> somaticSnvScoringModel;
    std::unique_ptr<VariantScoringModelServer> somaticIndelScoringModel;

private:
    std::unique_ptr<somatic_snv_caller_strand_grid> _sscaller_strand_grid;
    std::unique_ptr<somatic_indel_caller_grid> _sicaller_grid;
};

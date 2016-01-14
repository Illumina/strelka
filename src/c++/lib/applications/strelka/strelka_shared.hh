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
#include "calibration/VariantScoringModel.hh"
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
    double minimumEVS = 2.35;

    // TODO: STARKA-296 enable when ready
    bool is_use_indel_empirical_scoring = false;
};



struct strelka_options : public starling_base_options
{
    typedef starling_base_options base_t;

    strelka_options()
    {
        // turn on empirical scoring for strelka only:
        is_compute_somatic_scoring_metrics = true;
    }

    bool is_tumor_realigned_read() const
    {
        return (! tumor_realigned_read_filename.empty());
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

    std::string tumor_bam_filename;

    std::string tumor_realigned_read_filename;

    double somatic_snv_rate = 0.000001;
    std::string somatic_snv_filename;

    double somatic_indel_rate = 0.000001;
    std::string somatic_indel_filename;

    double shared_site_error_rate = 0.000005;
    double shared_site_error_strand_bias_fraction = 0.5;
    double site_somatic_normal_noise_rate = 0;
    bool is_site_somatic_normal_noise_rate = false;

    double shared_indel_error_factor = 1.65;
    double shared_indel_error_strand_bias_fraction = 0.1;
    double indel_somatic_normal_noise_rate = 0;
    bool is_indel_somatic_normal_noise_rate = false;

    // We provide a lower flank requirement for normal sample reads
    // during somatic variant calling, to ensure that all evidence
    // potentially used against a somatic call in the normal is
    // available:
    int normal_sample_min_read_bp_flank = 1;

    // Lower the candidate indel threshold on the tumor side to
    // increase sensitivity in case of low purity:
    //
    bool is_tumor_sample_min_candidate_indel_reads = false;
    bool is_tumor_sample_min_small_candidate_indel_read_frac = false;
    int tumor_sample_min_candidate_indel_reads = 2;
    double tumor_sample_min_small_candidate_indel_read_frac = 0.02;

    std::string somatic_callable_filename;

    // positions/indels in vcf are used to estimate low-frequency sequencing noise:
    std::vector<std::string> noise_vcf;

    somatic_filter_options sfilter;
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

    double max_chrom_depth = 0;
    double expected_chrom_depth = 0;
    cdmap_t chrom_depth;
    unsigned indelRegionStage = 0;
};



struct somatic_snv_caller_strand_grid;
struct somatic_indel_caller_grid;


// data deterministically derived from the input options, or read in from model files, etc.
//
struct strelka_deriv_options : public starling_base_deriv_options
{
    typedef starling_base_deriv_options base_t;

    strelka_deriv_options(
        const strelka_options& opt,
        const reference_contig_segment& ref);

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

    std::unique_ptr<VariantScoringModel> somaticSnvScoringModel;
    std::unique_ptr<VariantScoringModel> somaticIndelScoringModel;

private:
    std::unique_ptr<somatic_snv_caller_strand_grid> _sscaller_strand_grid;
    std::unique_ptr<somatic_indel_caller_grid> _sicaller_grid;
};

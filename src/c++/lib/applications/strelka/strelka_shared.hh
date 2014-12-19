// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Starka
// Copyright (c) 2009-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

///
/// \author Chris Saunders
///

#pragma once

#include "blt_util/chrom_depth_map.hh"
#include "starling_common/starling_shared.hh"


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

    unsigned indelMaxRefRepeat = 8;
    unsigned indelMaxIntHpolLength = 14;
    double indelMaxWindowFilteredBasecallFrac = 0.3;
    int sindelQuality_LowerBound = 30;

    unsigned indelRegionFlankSize = 50;
    double minimumQscore = 3.0;
};



struct strelka_options : public starling_options
{
    typedef starling_options base_t;

    bool is_tumor_bindel_diploid() const
    {
        return (! tumor_bindel_diploid_filename.empty());
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

    // report whether any type of indel-caller is running (including
    // checks from child class options):
    bool
    is_call_indels() const override
    {
        return (is_somatic_indel() || base_t::is_call_indels());
    }

    std::string tumor_bam_filename;

    std::string tumor_bindel_diploid_filename;
    std::string tumor_realigned_read_filename;

    double somatic_snv_rate = 0.000001;
    std::string somatic_snv_filename;

    double somatic_indel_rate = 0.000001;
    std::string somatic_indel_filename;

    double shared_site_error_rate = 0.000005;
    double shared_site_error_strand_bias_fraction = 0.5;
    double site_somatic_normal_noise_rate = 0;
    bool is_site_somatic_normal_noise_rate = false;

    double shared_indel_error_rate = 0.0000001;
    double shared_indel_error_factor = 1.2;
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

    double max_depth = 0;
    cdmap_t chrom_depth;
    unsigned indelRegionStage = 0;
};



struct somatic_snv_caller_strand_grid;
struct somatic_indel_caller_grid;


// data deterministically derived from the input options:
//
struct strelka_deriv_options : public starling_deriv_options
{
    typedef starling_deriv_options base_t;

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

private:
    std::unique_ptr<somatic_snv_caller_strand_grid> _sscaller_strand_grid;
    std::unique_ptr<somatic_indel_caller_grid> _sicaller_grid;
};

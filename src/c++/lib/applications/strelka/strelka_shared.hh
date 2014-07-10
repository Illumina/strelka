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

#include "starling_common/chrom_depth_map.hh"
#include "starling_common/starling_shared.hh"


/// variant call filtration options used only for somatic snvs and indels
///
/// these are isolated from starling gvcf options to avoid conflicts
///
struct somatic_filter_options
{
    somatic_filter_options()
        : is_skip_header(false)
        , max_depth_factor(3.)
        , snv_max_filtered_basecall_frac(0.4)
        , snv_max_spanning_deletion_frac(0.75)
        , snv_min_qss_ref(15)
        , indelMaxRefRepeat(8)
        , indelMaxIntHpolLength(14)
        , indelMaxWindowFilteredBasecallFrac(0.3)
        , sindelQuality_LowerBound(30)
        , indelRegionFlankSize(50)
    {}

    bool
    is_depth_filter() const
    {
        return (! chrom_depth_file.empty());
    }

    std::string chrom_depth_file;
    bool is_skip_header;
    double max_depth_factor;
    double snv_max_filtered_basecall_frac;
    double snv_max_spanning_deletion_frac;
    int snv_min_qss_ref;

    unsigned indelMaxRefRepeat;
    unsigned indelMaxIntHpolLength;
    double indelMaxWindowFilteredBasecallFrac;
    int sindelQuality_LowerBound;

    unsigned indelRegionFlankSize;
};



struct strelka_options : public starling_options
{
    typedef starling_options base_t;

    strelka_options()
        : somatic_snv_rate  (0.000001)
        , somatic_indel_rate(0.000001)
        , shared_site_error_rate(0.000005)
        , shared_site_error_strand_bias_fraction(0.5)
        , site_somatic_normal_noise_rate(0)
        , is_site_somatic_normal_noise_rate(false)
        , shared_indel_error_rate(0.0000001)
        , shared_indel_error_strand_bias_fraction(0.1)
        , indel_somatic_normal_noise_rate(0)
        , is_indel_somatic_normal_noise_rate(false)
        , normal_sample_min_read_bp_flank(1)
        , is_tumor_sample_min_candidate_indel_reads(false)
        , is_tumor_sample_min_small_candidate_indel_read_frac(false)
        , tumor_sample_min_candidate_indel_reads(2)
        , tumor_sample_min_small_candidate_indel_read_frac(0.02)
    {}

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
    virtual
    bool
    is_call_indels() const
    {
        return (is_somatic_indel() || base_t::is_call_indels());
    }

    std::string tumor_bam_filename;

    std::string tumor_bindel_diploid_filename;
    std::string tumor_realigned_read_filename;

    double somatic_snv_rate;
    std::string somatic_snv_filename;

    double somatic_indel_rate;
    std::string somatic_indel_filename;

    double shared_site_error_rate;
    double shared_site_error_strand_bias_fraction;
    double site_somatic_normal_noise_rate;
    bool is_site_somatic_normal_noise_rate;

    double shared_indel_error_rate;
    double shared_indel_error_strand_bias_fraction;
    double indel_somatic_normal_noise_rate;
    bool is_indel_somatic_normal_noise_rate;

    // We provide a lower flank requirement for normal sample reads
    // during somatic variant calling, to ensure that all evidence
    // potentially used against a somatic call in the normal is
    // available:
    int normal_sample_min_read_bp_flank;

    // Lower the candidate indel threshold on the tumor side to
    // increase sensitivity in case of low purity:
    //
    bool is_tumor_sample_min_candidate_indel_reads;
    bool is_tumor_sample_min_small_candidate_indel_read_frac;
    int tumor_sample_min_candidate_indel_reads;
    double tumor_sample_min_small_candidate_indel_read_frac;

    std::string somatic_callable_filename;

    somatic_filter_options sfilter;
};



/// somatic filter options computed after user input is finished:
///
struct somatic_filter_deriv_options
{
    somatic_filter_deriv_options()
        : max_depth(0.),
          indelRegionStage(0)
    {}

    bool
    is_max_depth() const
    {
        return (! chrom_depth.empty());
    }

    double max_depth;
    cdmap_t chrom_depth;
    unsigned indelRegionStage;
};



//struct somatic_snv_caller;
//struct somatic_snv_caller_grid;
struct somatic_snv_caller_strand_grid;
//struct somatic_indel_caller;
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

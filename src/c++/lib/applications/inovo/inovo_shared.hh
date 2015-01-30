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
#include "starling_common/starling_base_shared.hh"


/// variant call filtration options used only for denovo snvs and indels
///
/// these are isolated from starling gvcf options to avoid conflicts
///
struct denovo_filter_options
{
    bool
    is_depth_filter() const
    {
        return (! chrom_depth_file.empty());
    }

    /// this means *proband* depth
    std::string chrom_depth_file;
    bool is_skip_header = false;
    double max_depth_factor = 3.;
#if 0
    double snv_max_filtered_basecall_frac = 0.4;
    double snv_max_spanning_deletion_frac = 0.75;
    int snv_min_qss_ref = 15;

    unsigned indelMaxRefRepeat = 8;
    unsigned indelMaxIntHpolLength = 14;
    double indelMaxWindowFilteredBasecallFrac = 0.3;
    int sindelQuality_LowerBound = 30;

    unsigned indelRegionFlankSize = 50;
    double minimumQscore = 3.0;
#endif
};


namespace INOVO_SAMPLETYPE
{
    enum index_t
    {
        PROBAND,
        PARENT,
        SIBLING,
        SIZE
    };
}

namespace INOVO_GENDER
{
    enum index_t
    {
        UNKNOWN,
        MALE,
        FEMALE,
        SIZE
    };
}

/// tracks all sample information provided by the user
struct SampleInfo
{
    /// the id is only used to tell samples apart, this allows for future expansion beyond a 1-1 mapping of samples to bam files
    unsigned id = 0;

    /// relationship of sample to proband there's no use for an unknown value here:
    INOVO_SAMPLETYPE::index_t stype = INOVO_SAMPLETYPE::PROBAND;

    /// gender of sample, this is provided for future expansions but will be ignored for POC
    INOVO_GENDER::index_t gtype = INOVO_GENDER::UNKNOWN;
};



struct inovo_options : public starling_base_options
{
    typedef starling_base_options base_t;

    std::vector<std::string> alignmentFilename;

    std::vector<SampleInfo> alignmentSampleInfo;

#if 0
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

    double shared_indel_error_factor = 1.4;
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
#endif
    denovo_filter_options dfilter;
};



/// somatic filter options computed after user input is finished:
///
struct denovo_filter_deriv_options
{
    bool
    is_max_depth() const
    {
        return (! chrom_depth.empty());
    }

    double max_depth = 0;
    cdmap_t chrom_depth;
};



// data deterministically derived from the input options:
//
struct inovo_deriv_options : public starling_base_deriv_options
{
    typedef starling_base_deriv_options base_t;

    inovo_deriv_options(
        const inovo_options& opt,
        const reference_contig_segment& ref);

    ~inovo_deriv_options();

    denovo_filter_deriv_options dfilter;
};

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

#include "DenovoAlignmentFileOptions.hh"
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

    double dindel_qual_lowerbound = 30;
    double dsnv_qual_lowerbound = 30;
    double snv_max_filtered_basecall_frac = 0.35;
    double snv_max_spanning_deletion_frac = 0.75;
    int snv_min_qss_ref = 15;

    unsigned indelMaxRefRepeat = 14;
    unsigned indelMaxIntHpolLength = 14;
    double indelMaxWindowFilteredBasecallFrac = 0.3;
    unsigned sindelQuality_LowerBound = 5;

    unsigned indelRegionFlankSize = 50;
    double minimumQscore = 3.0;
};


struct pedicure_options : public starling_base_options
{
    typedef starling_base_options base_t;

    pedicure_options()
    {}

    bool
    is_denovo_callable() const
    {
        return (! denovo_callable_filename.empty());
    }


    DenovoAlignmentFileOptions alignFileOpt;

    /// variant vcf output file:
    std::string denovo_filename;

    /// callable bed track:
    std::string denovo_callable_filename;

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

    // report whether any type of indel-caller is running (including
    // checks from child class options):
    bool
    is_call_indels() const override
    {
        return (is_somatic_indel() || base_t::is_call_indels());
    }

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
    int tumor_sample_min_candidate_indel_reads = 2;
    // positions/indels in vcf are used to estimate low-frequency sequencing noise:
    std::vector<std::string> noise_vcf;
#endif
    denovo_filter_options dfilter;
};



/// denovo filter options computed after user input is finished:
///
struct denovo_filter_deriv_options
{
    bool
    is_max_depth() const
    {
        return (! chrom_depth.empty());
    }

    cdmap_t chrom_depth;
};



/// data deterministically derived from the input options
struct pedicure_deriv_options : public starling_base_deriv_options
{
    typedef starling_base_deriv_options base_t;

    explicit
    pedicure_deriv_options(const pedicure_options& opt);

    ~pedicure_deriv_options();

    denovo_filter_deriv_options dfilter;
};

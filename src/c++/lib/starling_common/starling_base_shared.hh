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

#include "starling_pos_processor_base_stages.hh"

#include "blt_common/blt_shared.hh"
#include "blt_util/PrettyFloat.hh"
#include "blt_util/reference_contig_segment.hh"
#include "starling_common/min_count_binom_gte_cache.hh"
#include "starling_common/starling_align_limit.hh"
#include "starling_common/Tier2Options.hh"

#include <cmath>


// starling max read size increases to the highest observed read size,
// up to the following practical limit:
//
enum { STARLING_MAX_READ_SIZE = 25000 };



struct starling_base_options : public blt_options
{
    typedef blt_options base_t;

    starling_base_options()
//      : readConfidentSupportThreshold(readConfidentSupportThresholdDefault)
    {}

    virtual
    bool
    is_all_sites() const
    {
        return false;
    }

    bool
    is_write_candidate_indels() const
    {
        return (! candidate_indel_filename.empty());
    }

    // parameters inherited from varling caller:
    //
    double bindel_diploid_theta = 0.0001;

    // use this theta in long homopolymers (germline only)
    double indelHighRepeatTheta = 0.01;
    // definition of "long homopolymer" above
    unsigned indelHighRepeatCount = 16;

    uint32_t user_genome_size = 0; // genome size specified by user for the indel calling model -- actual value used is in deriv_options.
    bool is_user_genome_size = false;

    // parameter to enable/disable short haplotype calling
    // it's always false for now
    bool is_short_haplotype_calling_enabled = true;

    // to contribute to a breakpoint likelihood, a read must have at least
    // this many bases on each side of the breakpoint:
    //
    // This is the default used in all samples unless an override is provided for the sample.
    //
    int default_min_read_bp_flank = 5;

    // starling parameters:
    //

    // should reads falling below the snp-caller's mapping criteria be
    // realigned? (note this only makes sense if writing out realigned
    // reads
    bool is_realign_submapped_reads = false;

    // maximum indel size which can be represented by starling
    // (formerly a static value)
    unsigned max_indel_size = 150;

    // Do we test indel observation counts to determine if these are significant enough
    // to create an indel candidate? This should be true for any normal variant caller,
    // it is turned off for specialized indel noise estimation routines
    bool is_candidate_indel_signal_test = true;

    // Observed indels are promoted to candidate indels based on a one-sided
    // binomial exact test, which incorporates expected per-read indel error rate,
    // total coverage, and observed indel coverage.
    //
    // this sets the p-value threshold for determining indel candidacy, a lower value
    // means that relatively more observations are required to create any indel candidate
    const double indel_candidate_signal_test_alpha = 1e-9;

    int max_read_indel_toggle = 5; // if a read samples more than max indel changes, we skip realignment
    double max_candidate_indel_density = 0.15; // max number of candidate indels per read base, if exceeded search is curtailed to toggle depth=1

    // max estimated read depth for indels to reach candidacy for realignment and indel calling. Non-positive value disables this
    int max_candidate_indel_depth = -1;

    // depth factor above filtration filter cutoff where
    // indel candidacy is disallowed:
    double max_candidate_indel_depth_factor = 3.;

    // min length of the 'inserted' portion of an open-ended breakpoint:
    unsigned min_candidate_indel_open_length = 20;

    // the maximum number of candidate re-alignments for each read:
    unsigned max_realignment_candidates = 5000;

    // clip the section of a read which aligns equally well to two or
    // more paths before pileup or realigned read output
    bool is_clip_ambiguous_path = true;

    bool is_realigned_read_file = false;

    // this option imposes a consistency criteria on alignments with
    // nearly equal score to favor certain alignments even if they do
    // not have the optimal score.
    //
    // when using smoothed alignments, we realign to the prefered
    // alignment that is not more than smoothed_lnp_range from the
    // highest scoring alignment.
    //
    bool is_smoothed_alignments = true;
    double smoothed_lnp_range = std::log(10.);

    // filter reads where both reads of pair have SE score 0, temp fix
    // for internal analysis:
    bool is_filter_unanchored = false;

    // some newer applications will not use the "bam_filename" input:
    bool is_bam_filename_used = true;

    std::string realigned_read_filename;
    std::string bam_filename; // BAM/CRAM input file
    std::string bam_seq_name;

    std::string candidate_indel_filename;
    bool is_write_candidate_indels_only = false;

    double indel_nonsite_match_prob = 0.25;

    // Assume a ploidy-based prior (0%, 50%, 100% or 0% 100% for haploid, diploid
    // If false, a continuous model is used
    bool is_ploidy_prior = true;

    // the assumed noise (uniform) in basecalls
    double noise_floor = -1;

    // the minimum allele frequency to call a heterozygous genotype when not assuming ploidy
    double min_het_vf = 0.01;

    // vcfs can be input to specify candidate indels:
    std::vector<std::string> input_candidate_indel_vcf;

    // positions/indels in vcf must be written in output:
    std::vector<std::string> force_output_vcf;

    // Indicates that an upstream oligo is present on reads, which can be used to increase confidence for indels near the edge of the read
    unsigned upstream_oligo_size = 0;

    /// file specifying haploid and deleted regions as 1/0 in bed col 4
    std::string ploidy_region_bedfile;

    /// Stores runtime stats
    std::string segmentStatsFilename;

    bool
    isMaxBufferedReads() const
    {
        return (maxBufferedReads != 0);
    }

    /// maximum number of reads which can be buffered for any one sample
    ///
    /// set to zero to disable limit
    unsigned maxBufferedReads = 100000;

    bool isBasecallQualAdjustedForMapq = true;

    Tier2Options tier2;

    // indel error options
    std::string indel_error_models_filename;
    std::string indel_error_model_name = "logLinear";

    // temporary indel erorr rate hack applied to germline only
    bool isIndelErrorRateFactor = false;
    double indelErrorRateFactor = 1.;

    // when P(read | allele) is >= this value the read counts as "supporting"
    // the given allele
    // WARNING: this value impacts several count-based EVS metrics
    PrettyFloat<double> readConfidentSupportThreshold = PrettyFloat<double>("0.9");
};



// allow for sample-specific parameter values:
//
struct starling_sample_options
{
    explicit
    starling_sample_options(
        const starling_base_options& opt)
        : min_read_bp_flank(opt.default_min_read_bp_flank)
    {}

    int min_read_bp_flank;
};



struct IndelErrorModel;
struct indel_digt_caller;
struct GenotypePriorSet;


// data deterministically derived from the input options:
//
struct starling_base_deriv_options : public blt_deriv_options
{
    typedef blt_deriv_options base_t;

    starling_base_deriv_options(
        const starling_base_options& opt,
        const reference_contig_segment& ref);

    ~starling_base_deriv_options();

    const indel_digt_caller&
    incaller() const
    {
        return *(_incaller.get());
    }

    double
    get_nonsite_path_lnp(const bool is_tier2_pass,
                         const uint16_t nsite) const
    {
        const double nonsite_match_lnp(is_tier2_pass ?
                                       tier2_indel_nonsite_match_lnp :
                                       indel_nonsite_match_lnp );
        return nonsite_match_lnp*nsite;
    }

    const std::vector<unsigned>&
    getPostCallStage() const
    {
        return _postCallStage;
    }

    const IndelErrorModel&
    getIndelErrorModel() const
    {
        return *_indelErrorModel;
    }

    const GenotypePriorSet&
    getIndelGenotypePriors() const
    {
        return *_indelGenotypePriors;
    }

protected:
    unsigned
    addPostCallStage(
        const unsigned size)
    {
        _postCallStage.push_back(size);
        return STAGE::SIZE+_postCallStage.size()-1;
    }

public:
    double indel_nonsite_match_lnp;
    double tier2_indel_nonsite_match_lnp;

    double site_lnprior;
    double nonsite_lnprior;

    starling_align_limit sal;

    unsigned variant_window_first_stage;
    unsigned variant_window_last_stage;

    const min_count_binom_gte_cache countCache;

    // cache the log of options->indelErrroRateFactor
    const double logIndelErrorRateFactor;

private:
    std::unique_ptr<IndelErrorModel> _indelErrorModel;
    std::unique_ptr<indel_digt_caller> _incaller; // object to precalculate bindel_diploid priors..
    std::unique_ptr<GenotypePriorSet> _indelGenotypePriors;

    std::vector<unsigned> _postCallStage;
};



struct starling_read_counts : public blt_read_counts
{
    void
    report(std::ostream& os) const;

    unsigned normal_indel_used = 0;
    unsigned normal_indel_intersect = 0;
};

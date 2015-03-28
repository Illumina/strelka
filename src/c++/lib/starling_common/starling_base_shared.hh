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

#include "starling_pos_processor_base_stages.hh"

#include "blt_common/blt_shared.hh"
#include "blt_util/reference_contig_segment.hh"
#include "starling_common/starling_align_limit.hh"
#include "starling_common/Tier2Options.hh"

#include <cmath>


// starling max read size increases to the highest observed read size,
// up to the following practical limit:
//
enum { STARLING_MAX_READ_SIZE = 25000 };



struct avg_window_data
{
    bool
    operator<(const avg_window_data& rhs) const
    {
        return (flank_size < rhs.flank_size);
    }

    unsigned flank_size;
    std::string filename;
};

std::ostream&
operator<<(std::ostream& os, const avg_window_data& awd);



struct starling_base_options : public blt_options
{
    typedef blt_options base_t;

    starling_base_options() {}

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

    unsigned htype_buffer_segment() const
    {
        return max_indel_size;
    }

    // parameters inherited from varling caller:
    //
    double bindel_diploid_theta = 0.0001;
    double bindel_diploid_het_bias = 0;
    bool is_bindel_diploid_het_bias = false;
    bool is_test_indels = false;
    uint32_t user_genome_size = 0; // genome size specified by user for the indel calling model -- actual value used is in deriv_options.
    bool is_user_genome_size = false;

    /// to contribute to a breakpoint likelihood, a read must have at least
    /// this many bases on each side of the breakpoint:
    ///
    /// This is the default used in all samples unless an override is provided for the sample.
    ///
    int default_min_read_bp_flank = 6;

    // starling parameters:
    //

    // should reads falling below the snp-caller's mapping criteria be
    // realigned? (note this only makes sense if writing out realigned
    // reads
    bool is_realign_submapped_reads = false;

    // maximum indel size which can be represented by starling
    // (formerly a static value)
    unsigned max_indel_size = 150;

    // indel cannot become candidate unless at least min reads which
    // meet mapping threshold support it
    //
    // this is the default used for all samples until overridden
    int default_min_candidate_indel_reads = 3;

    // indel cannot become candidate unless at least frac of reads
    // which meet mapping thresholds support it (num is reads
    // supporting indel/den in aligner reads aligning to adjacent
    // position).
    double min_candidate_indel_read_frac = 0.02;

    // indels this size or lower have additional 'small indel'
    // candidacy criteria
    int max_small_candidate_indel_size = 4;

    // same as for min_candidate_indel_read_frac, but for small indels
    //
    // this is the default used for all samples until overridden
    double default_min_small_candidate_indel_read_frac = 0.1;

    // indels in homopolymers use a one-sided binomial exact test to
    // determine whether they are eligible for candidacy, based on the
    // expected per-read error rate, total coverage, and indel coverage
    // this sets the p-value threshold for determining homopolymer candidacy
    double tumor_min_hpol_pval = std::pow(10,-9);

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

    std::vector<avg_window_data> variant_windows;

    // Test if an indel is not in the two most likely indel alleles
    // among a conflicting set.  if not remove from variant
    // calling. This option is for single sample calling only.
    //
    bool is_noise_indel_filter = false;

    // vcfs can be input to specify candidate indels:
    std::vector<std::string> input_candidate_indel_vcf;

    // positions/indels in vcf must be written in output:
    std::vector<std::string> force_output_vcf;

    // Internal development option - not for production use:
    bool is_baby_elephant = false;

    // if not using grouper, you can optionally turn off the restriction that each qname occurs once in the bam:
    bool is_ignore_read_names = false;

    // Indicates that an upstream oligo is present on reads, which can be used to increase confidence for indels near the edge of the read
    unsigned upstream_oligo_size = 0;

    // turn on haplotype based calling:
    bool is_htype_calling = false;

    // number of allowed haplotypes:
    unsigned hytpe_count = 2;

    // multiple of max-indel-size used for haplotyping:
    unsigned htype_call_segment = 1000;

    // if true, treat all soft-clipped segments on the egdes of reads as realignable
    bool is_remap_input_softclip = true;

    /// file specifying haploid and deleted regions as 1/0 in bed col 4
    std::string ploidy_region_bedfile;

    bool
    isMaxBufferedReads() const
    {
        return (maxBufferedReads != 0);
    }

    /// maximum number of reads which can be buffered for any one sample
    ///
    /// set to zero to disable limit
    unsigned maxBufferedReads = 100000;

    Tier2Options tier2;
};



// allow for sample-specific parameter values:
//
struct starling_sample_options
{
    starling_sample_options(
        const starling_base_options& opt)
        : min_read_bp_flank(opt.default_min_read_bp_flank)
        , min_candidate_indel_reads(opt.default_min_candidate_indel_reads)
        , min_small_candidate_indel_read_frac(opt.default_min_small_candidate_indel_read_frac)
    {}

    int min_read_bp_flank;
    int min_candidate_indel_reads;
    double min_small_candidate_indel_read_frac;
};



struct indel_digt_caller;



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
    getPostCalLStage() const
    {
        return _postCallStage;
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

    std::string bam_header_data; // the full bam header, read in from bam file. Used for setting the sample name in

    unsigned variant_window_first_stage;
    unsigned variant_window_last_stage;

private:
    std::unique_ptr<indel_digt_caller> _incaller; // object to precalculate bindel_diploid priors..

    std::vector<unsigned> _postCallStage;
};



struct starling_read_counts : public blt_read_counts
{
    void
    report(std::ostream& os) const;

    unsigned normal_indel_used = 0;
    unsigned normal_indel_intersect = 0;
};

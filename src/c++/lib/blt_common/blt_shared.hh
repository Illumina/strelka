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

#include "blt_util/blt_types.hh"
#include "blt_util/pos_range.hh"
#include "blt_util/seq_util.hh"

#include <memory>
#include <vector>

extern const char STDIN_FILENAME[];

extern const unsigned MAX_FLANK_SIZE;


namespace LOG_LEVEL
{
enum index_t
{
    DEFAULT = 0, // 0 = errors and low-frequency warnings
    ALLWARN = 1  // all other warnings
};
}



struct blt_options
{
    blt_options() {}
    virtual ~blt_options() {}

    bool
    is_ref_set() const
    {
        return (is_samtools_ref_set);
    }

    bool
    is_nonref_test() const
    {
        return (! nonref_test_filename.empty());
    }

    bool
    is_nonref_sites() const
    {
        return (! nonref_sites_filename.empty());
    }

    bool
    is_compute_germline_VQSRmetrics() const
    {
        return (is_report_germline_VQSRmetrics || (! calibration_model.empty()));
    }

    virtual
    bool
    is_bsnp_diploid() const { return false; }

    bool
    is_tier2() const
    {
        return
            (is_tier2_min_single_align_score ||
             is_tier2_min_paired_align_score ||
             is_tier2_single_align_score_rescue_mode ||
             is_tier2_mismatch_density_filter_count ||
             is_tier2_no_mismatch_density_filter ||
             is_tier2_no_filter_unanchored ||
             is_tier2_include_singleton ||
             is_tier2_include_anomalous);
    }

    double lsnp_alpha = 0;
    double bsnp_diploid_theta = 0.001;
    double bsnp_monoploid_theta = 0;
    int bsnp_nploid_ploidy = 0;
    double bsnp_nploid_snp_prob = 0;
    double bsnp_ssd_no_mismatch = 0;
    double bsnp_ssd_one_mismatch = 0;
    double bsnp_diploid_het_bias = 0;

    double adis_lrt_alpha = 0;
    double adis_table_alpha = 0;
    double adis_win_lrt_alpha = 0;
    unsigned adis_win_lrt_flank_size = 0;
    double acov_alpha = 0;
    bool is_lsnp = false;
    bool is_bsnp_monoploid = false;
    bool is_bsnp_nploid = false;
    bool is_bsnp_diploid_het_bias = false;
    bool is_adis_lrt = false;
    bool is_adis_table = false;
    bool is_adis_win_lrt = false;
    bool is_acov = false;

    int min_qscore = 17;
    int min_single_align_score = 10;
    int min_paired_align_score = 6;
    bool single_align_score_exclude_mode = false;
    bool single_align_score_rescue_mode = false;

    int tier2_min_single_align_score = 0;
    bool is_tier2_min_single_align_score = false;
    int tier2_min_paired_align_score = 0;
    bool is_tier2_min_paired_align_score = false;
    bool is_tier2_single_align_score_rescue_mode = false;

    int tier2_mismatch_density_filter_count = 0;
    bool is_tier2_mismatch_density_filter_count = false;

    bool is_tier2_no_mismatch_density_filter = false;
    bool is_tier2_no_filter_unanchored = false;
    bool is_tier2_include_singleton = false;
    bool is_tier2_include_anomalous = false;

    bool is_min_win_qscore = false;
    int min_win_qscore = 0;
    unsigned min_win_qscore_flank_size = 0;
    bool is_max_win_mismatch = false;
    unsigned max_win_mismatch = 0;
    unsigned max_win_mismatch_flank_size = 0;
    bool is_counts = false;
    bool is_print_evidence = false;
    bool is_print_all_site_evidence = false;
    pos_range user_report_range;   // requested report range
    bool is_read_sample = false;
    double read_sample_rate = 0;

    bool is_samtools_ref_set = false;
    std::string samtools_ref_seq_file;

    bool is_filter_anom_calls = false;
    bool is_include_singleton = false;
    bool is_include_anomalous = false;

    std::string counts_filename;

    bool is_clobber = true;
    bool is_report_range_ref = false;
    bool is_print_all_poly_gt = false; // print the posterior probabilities for all genotypes
    bool is_print_used_allele_counts = false; // print allele counts as in CASAVA 1.7 output
    int used_allele_count_min_qscore = 0; // print the above with a qscore cutoff...
    double max_basecall_filter_fraction = 1.; // if more than this fraction of basecalls are filtered out, than filter the snp

    int max_vexp_iterations = 0;
    bool is_min_vexp = false;
    double min_vexp = 0;

    LOG_LEVEL::index_t verbosity = LOG_LEVEL::DEFAULT;

    bool is_write_variable_metadata = true;

    std::string cmdline;

    // constants for het-bias model:
    //
    // Humans will often pick exact multples of the max_ratio increment,
    // which are also the least efficient points in terms of increment
    // size -- fudge removes this trend from the computation:
    //
    static constexpr double het_bias_inc_fudge = 0.0001;
    const double het_bias_max_ratio_inc = 0.05 + het_bias_inc_fudge;

    double nonref_variant_rate = 0.000001;
    double min_nonref_freq = 0;
    double nonref_site_error_rate = 0.0001;
    double nonref_site_error_decay_freq = 0.01;
    std::string nonref_test_filename;
    std::string nonref_sites_filename;

    bool is_eland_compat = false;

    bool is_max_input_depth = false;
    unsigned max_input_depth = 0;

    bool is_compute_hapscore = false;
    bool is_report_germline_VQSRmetrics = false;
    bool is_compute_calibration_features = false;// For development only, out all features needed im

    bool is_compute_somatic_VQSRmetrics = false;

    // Which calibration model should we use?
    // leave blank to select default rule-based metric option
    std::string calibration_model;

    // Apply codon phasing:
    bool do_codon_phasing = false;

    // Size of the window we are phasing in, default is codon range (=3)
    int phasing_window = 3;

    // Apply assembly to qualifying regions
    bool do_assemble                        = false;
    unsigned assemble_aggresiveness         = 1;
    std::string assembly_regions_filename;

    //multiplier for ref error rate to use; 1 would be expected to be correct, but higher values counteract a bias away from homozygous indels (undercalling)
    double indel_ref_error_factor = 1.;

    std::string report_filename;
    std::string calibration_models_filename;
    std::string indel_scoring_models;   // file containing all indel scoring models
    std::string indel_error_model;      // which baseline prior should be used for candidate indel genotyping
};



struct pprob_digt_caller;


// data deterministically derived from the user input options:
//
struct blt_deriv_options
{
    /// @param ref_end this is either the full reference contig size,
    /// or the end position of the acquired reference segment if
    /// -report-range-end was used
    ///
    blt_deriv_options(
        const blt_options& opt,
        const pos_t ref_end);

    ~blt_deriv_options();

    pos_range report_range;
    pos_range report_range_limit;   //  maximum report range

    const pprob_digt_caller&
    pdcaller() const
    {
        return *(_pdcaller.get());
    }

private:
    std::unique_ptr<pprob_digt_caller> _pdcaller; // object to precalculate bsnp_diploid priors..
};



struct blt_read_counts
{
    void
    report(std::ostream& os) const;

    unsigned subsample_filter = 0;
    unsigned primary_filter = 0;
    unsigned duplicate = 0;
    unsigned unmapped = 0;
    unsigned secondary = 0;
    unsigned supplement = 0;
    unsigned unanchored = 0;
    unsigned large_ref_deletion = 0;
    unsigned align_score_filter = 0;

    // Floating means the read is indicated as mapped, but has no "M"
    // in the cigar string. Typically inside of an insertion.
    unsigned floating = 0;

    // if optional setting is given to filter out reads once a certain depth
    // is exceeded, the number of reads filtered are enumerated here:
    unsigned max_depth = 0;
    unsigned used = 0;
};



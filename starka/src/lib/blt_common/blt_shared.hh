// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2012 Illumina, Inc.
//
// This software is covered by the "Illumina Genome Analyzer Software
// License Agreement" and the "Illumina Source Code License Agreement",
// and certain third party copyright/licenses, and any user of this
// source file is bound by the terms therein (see accompanying files
// Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
// Illumina_Source_Code_License_Agreement.pdf and third party
// copyright/license notices).
//
//

/// \file

/// \author Chris Saunders
///
#ifndef __BLT_SHARED_HH
#define __BLT_SHARED_HH

#include "blt_util/blt_types.hh"
#include "blt_util/pos_range.hh"
#include "blt_util/seq_util.hh"

#include <memory>
#include <vector>

extern const char STDIN_FILENAME[];

extern const unsigned MAX_FLANK_SIZE;

enum { BLT_MAX_READ_SIZE = 250 };

enum { MAX_READ_REF_DELETION_SIZE = 150 };


namespace LOG_LEVEL {
    enum index_t {
        DEFAULT = 0, // 0 = errors and low-frequency warnings
        ALLWARN = 1  // all other warnings
    };
}


struct gvcf_options {

    gvcf_options()
        : is_skip_header(false)
        , is_max_depth_factor(true)
        , max_depth_factor(3.)
        , is_min_gqx(true)
        , min_gqx(30.)
        , is_max_base_filt(true)
        , max_base_filt(.3)
        , is_max_snv_sb(true)
        , max_snv_sb(10)
        , is_max_snv_hpol(true)
        , max_snv_hpol(6)
        , is_max_ref_rep(true)
        , max_ref_rep(8)
        , block_label("BLOCKAVG_min30p3a")
        , block_frac_tol(.3)
        , block_abs_tol(3)
        , block_max_nonref(.2)
    {}

    // admin/other:
    std::string out_file;
    std::string chrom_depth_file;
    bool is_skip_header;

    // filters:
    bool is_max_depth_factor;
    double max_depth_factor;
    bool is_min_gqx;
    double min_gqx;
    bool is_max_base_filt;
    double max_base_filt;
    bool is_max_snv_sb;
    double max_snv_sb;
    bool is_max_snv_hpol;
    unsigned max_snv_hpol;
    bool is_max_ref_rep;
    unsigned max_ref_rep;

    // blocking scheme:
    std::string block_label;
    double block_frac_tol;
    int block_abs_tol;

    double block_max_nonref; // what percentage of non-ref bases can a site have and still be included in a non-variant block
};
    



struct blt_options {

    blt_options()
        : lsnp_alpha(0),
          bsnp_diploid_theta(0.001),
          bsnp_monoploid_theta(0),
          bsnp_nploid_ploidy(0),
          bsnp_nploid_snp_prob(0),
          bsnp_ssd_no_mismatch(0),
          bsnp_ssd_one_mismatch(0),
          bsnp_diploid_het_bias(0),
          adis_lrt_alpha(0),
          adis_table_alpha(0),
          adis_win_lrt_alpha(0),
          adis_win_lrt_flank_size(0),
          acov_alpha(0),
          is_lsnp(false),
          is_bsnp_monoploid(false),
          is_bsnp_nploid(false),
          is_bsnp_diploid_file(false),
          is_bsnp_diploid_allele_file(false),
          is_bsnp_diploid_het_bias(false),
          is_adis_lrt(false),
          is_adis_table(false),
          is_adis_win_lrt(false),
          is_acov(false),
          min_qscore(17),
          min_single_align_score(10),
          min_paired_align_score(6),
          single_align_score_exclude_mode(false),
          single_align_score_rescue_mode(false),

          tier2_min_single_align_score(0),
          is_tier2_min_single_align_score(false),
          tier2_min_paired_align_score(0),
          is_tier2_min_paired_align_score(false),
          is_tier2_single_align_score_rescue_mode(false),

          tier2_mismatch_density_filter_count(0),
          is_tier2_mismatch_density_filter_count(false),

          is_tier2_no_mismatch_density_filter(false),
          is_tier2_no_filter_unanchored(false),
          is_tier2_include_singleton(false),
          is_tier2_include_anomalous(false),

          is_min_win_qscore(false),
          min_win_qscore(0),
          min_win_qscore_flank_size(0),
          is_max_win_mismatch(false),
          max_win_mismatch(0),
          max_win_mismatch_flank_size(0),
          is_counts(false),
          is_print_evidence(false),
          is_print_all_site_evidence(false),
          is_read_sample(false),
          read_sample_rate(0),
          is_samtools_ref_set(false),
          is_filter_anom_calls(false),
          is_include_singleton(false),
          is_include_anomalous(false),
          is_clobber(false),
          is_report_range_ref(false),
          is_print_all_poly_gt(false),
          is_print_used_allele_counts(false),
          used_allele_count_min_qscore(0),
          max_basecall_filter_fraction(1.),
          max_vexp_iterations(0),
          is_min_vexp(false),
          min_vexp(0),
          verbosity(LOG_LEVEL::DEFAULT)
        , is_write_variable_metadata(true)
        , het_bias_inc_fudge(0.0001)
        , het_bias_max_ratio_inc(0.05+het_bias_inc_fudge)
        , nonref_variant_rate(0.000001)
        , min_nonref_freq(0)
        , nonref_site_error_rate(0.0001)
        , nonref_site_error_decay_freq(0.01)
        , is_eland_compat(false)
        , is_max_input_depth(false)
        , max_input_depth(0)
        , is_compute_hapscore(true)
    {}

    virtual ~blt_options() {}

    bool
    is_ref_set() const { return (is_samtools_ref_set); }

    bool
    is_nonref_test() const { return (! nonref_test_filename.empty()); }

    bool
    is_nonref_sites() const { return (! nonref_sites_filename.empty()); }

    // is the diploid snp model being used?
    bool
    is_bsnp_diploid() const {
        return (is_bsnp_diploid_file ||
                is_bsnp_diploid_allele_file ||
                is_gvcf_output());
    }

    bool
    is_all_sites() const {
        return (is_bsnp_diploid_allele_file || is_gvcf_output());
    }

    bool
    is_gvcf_output() const {
        return (! gvcf.out_file.empty());
    }

    bool
    is_tier2() const {
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

    double lsnp_alpha;
    double bsnp_diploid_theta;
    double bsnp_monoploid_theta;
    int bsnp_nploid_ploidy;
    double bsnp_nploid_snp_prob;
    double bsnp_ssd_no_mismatch;
    double bsnp_ssd_one_mismatch;
    double bsnp_diploid_het_bias;
    double adis_lrt_alpha;
    double adis_table_alpha;
    double adis_win_lrt_alpha;
    unsigned adis_win_lrt_flank_size;
    double acov_alpha;
    bool is_lsnp;
    bool is_bsnp_monoploid;
    bool is_bsnp_nploid;
    bool is_bsnp_diploid_file;
    bool is_bsnp_diploid_allele_file;
    bool is_bsnp_diploid_het_bias;
    bool is_adis_lrt;
    bool is_adis_table;
    bool is_adis_win_lrt;
    bool is_acov;
    int min_qscore;
    int min_single_align_score;
    int min_paired_align_score;
    bool single_align_score_exclude_mode;
    bool single_align_score_rescue_mode;

    int tier2_min_single_align_score;
    bool is_tier2_min_single_align_score;
    int tier2_min_paired_align_score;
    bool is_tier2_min_paired_align_score;
    bool is_tier2_single_align_score_rescue_mode;

    int tier2_mismatch_density_filter_count;
    bool is_tier2_mismatch_density_filter_count;

    bool is_tier2_no_mismatch_density_filter;
    bool is_tier2_no_filter_unanchored;
    bool is_tier2_include_singleton;
    bool is_tier2_include_anomalous;

    bool is_min_win_qscore;
    int min_win_qscore;
    unsigned min_win_qscore_flank_size;
    bool is_max_win_mismatch;
    unsigned max_win_mismatch;
    unsigned max_win_mismatch_flank_size;
    bool is_counts;
    bool is_print_evidence;
    bool is_print_all_site_evidence;
    pos_range user_report_range;   // requested report range
    bool is_read_sample;
    double read_sample_rate;

    bool is_samtools_ref_set;
    std::string samtools_ref_seq_file;

    bool is_filter_anom_calls;
    bool is_include_singleton;
    bool is_include_anomalous;

    std::string counts_filename;
    std::string bsnp_diploid_filename;
    std::string bsnp_diploid_allele_filename;

    bool is_clobber;
    bool is_report_range_ref;
    bool is_print_all_poly_gt; // print the posterior probabilities for all genotypes
    bool is_print_used_allele_counts; // print allele counts as in CASAVA 1.7 output
    int used_allele_count_min_qscore; // print the above with a qscore cutoff...
    double max_basecall_filter_fraction; // if more than this fraction of basecalls are filtered out, than filter the snp
    int max_vexp_iterations;
    bool is_min_vexp;
    double min_vexp;


    LOG_LEVEL::index_t verbosity;

    bool is_write_variable_metadata;

    std::string cmdline;

    // constants for het-bias model:
    //
    // Humans will often pick exact multples of the max_ratio increment,
    // which are also the least efficient points in terms of increment
    // size -- fudge removes this trend from the computation:
    //
    const double het_bias_inc_fudge;
    const double het_bias_max_ratio_inc;

    double nonref_variant_rate;
    double min_nonref_freq;
    double nonref_site_error_rate;
    double nonref_site_error_decay_freq;
    std::string nonref_test_filename;
    std::string nonref_sites_filename;

    bool is_eland_compat;

    bool is_max_input_depth;
    unsigned max_input_depth;

    bool is_compute_hapscore;

    std::string report_filename;

    gvcf_options gvcf;
};



struct pprob_digt_caller;


// data deterministically derived from the user input options:
//
struct blt_deriv_options {

    /// @param ref_end this is either the full reference contig size,
    /// or the end position of the acquired reference segment if
    /// -report-range-end was used
    ///
    blt_deriv_options(const blt_options& opt,
                      const pos_t ref_end);

    ~blt_deriv_options();

    pos_range report_range;
    pos_range report_range_limit;   //  maximum report range

    const pprob_digt_caller&
    pdcaller() const { return *(_pdcaller.get()); }

private:
    std::auto_ptr<pprob_digt_caller> _pdcaller; // object to precalculate bsnp_diploid priors..
};



struct blt_read_counts {

    blt_read_counts()
        : subsample_filter(0),
          primary_filter(0),
          duplicate(0),
          unmapped(0),
          secondary(0),
          unanchored(0),
          large_ref_deletion(0),
          align_score_filter(0),
          floating(0),
          max_depth(0),
          used(0) {}

    void
    report(std::ostream& os) const;

    unsigned subsample_filter;
    unsigned primary_filter;
    unsigned duplicate;
    unsigned unmapped;
    unsigned secondary;
    unsigned unanchored;
    unsigned large_ref_deletion;
    unsigned align_score_filter;

    // Floating means the read is indicated as mapped, but has no "M"
    // in the cigar string. Typically inside of an insertion.
    unsigned floating;

    // if optional setting is given to filter out reads once a certain depth
    // is exceeded, the number of reads filtered are enumerated here:
    unsigned max_depth;
    unsigned used;
};


#endif

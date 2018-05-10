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

#pragma once

#include "starling_pos_processor_base_stages.hh"

#include "blt_common/blt_shared.hh"

#include "blt_util/PrettyFloat.hh"
#include "blt_util/reference_contig_segment.hh"
#include "options/AlignmentFileOptions.hh"
#include "starling_common/min_count_binom_gte_cache.hh"
#include "starling_common/starling_align_limit.hh"
#include "starling_common/Tier2Options.hh"

#include <cmath>


/// This max read size is used as a QC check.
enum { STRELKA_MAX_READ_SIZE = 25000 };


typedef std::vector<std::string> regions_t;


struct starling_base_options : public blt_options
{
    typedef blt_options base_t;

    starling_base_options()
    {}

    void
    validate() const
    {
        base_t::validate();

        assert((min_het_vf > 0.) && (min_het_vf < 0.5));
        tier2.validate();
    }

    virtual
    bool
    is_all_sites() const
    {
        return false;
    }

    bool
    isWriteRealignedReads() const
    {
        return (not realignedReadFilenamePrefix.empty());
    }

    bool
    isUseCallRegions() const
    {
        return (not callRegionsBedFilename.empty());
    }

    /// \brief Provide access the input alignment file list
    virtual
    const AlignmentFileOptions&
    getAlignmentFileOptions() const = 0;

    unsigned
    getSampleCount() const
    {
        return getAlignmentFileOptions().alignmentFilenames.size();
    }

    std::string referenceFilename;

    // list of chromosome regions to be analyzed
    regions_t regions;

    //
    double bindel_diploid_theta = 0.0001;

    /// parameter to enable/disable short haplotype calling
    bool isHaplotypingEnabled = false;

    /// true if we run somatic calling and false otherwise
    bool isSomaticCallingMode = false;

    // to contribute to a breakpoint likelihood, a read must have at least
    // this many bases on each side of the breakpoint:
    //
    // This is the default used in all samples unless an override is provided for the sample.
    //
    int default_min_read_bp_flank = 5;

    // starling parameters:
    //

    /// If true, reads which fail the variant calling mapping thresholds are realigned using
    /// the same procedure as the variant calling reads
    bool is_realign_submapped_reads = false;

    /// \brief Maximum indel size.
    ///
    /// This is the maxumum size of an indel for both internal representation and external reporting.
    ///
    /// Indels detected above this size may still be reflected within the variant analysis process via the caller's
    /// open breakpoint mechanism.
    ///
    unsigned maxIndelSize = 49;

    // Do we test indel observation counts to determine if these are significant enough
    // to create an indel candidate? This should be true for any normal variant caller,
    // it is turned off for specialized indel noise estimation routines
    bool is_candidate_indel_signal_test = true;

    /// Observed indels are promoted to candidate indels based on a one-sided
    /// binomial exact test, which incorporates expected per-read indel error rate,
    /// total coverage, and observed indel coverage.
    ///
    /// this sets the p-value threshold for determining indel candidacy, a lower value
    /// means that relatively more observations are required to create any indel candidate
    const double indel_candidate_signal_test_alpha = 1e-9;

    int max_read_indel_toggle = 5; // if a read samples more than max indel changes, we skip realignment

    /// If there are more than this many candidate indels per base intersecting a read, then realignment
    /// is truncated to only allow individual indel toggles of the starting alignments for that read.
    ///
    /// Value must be greater than 0
    double max_candidate_indel_density = 0.15;

    // max estimated read depth for indels to reach candidacy for realignment and indel calling. If estimated
    // depth summed over all (non-tumor) samples exceeds this value in any sample, candidacy is disabled.
    // Non-positive value turns this filter off.
    int max_candidate_indel_depth = -1;

    // depth factor above filtration filter cutoff where
    // indel candidacy is disallowed:
    double max_candidate_indel_depth_factor = 3.;

    // min length of the 'inserted' portion of an open-ended breakpoint:
    unsigned min_candidate_indel_open_length = 20;

    // the maximum number of candidate re-alignments for each read:
    unsigned max_realignment_candidates = 5000;

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

    /// Path prefix for all realigned bam output files
    std::string realignedReadFilenamePrefix;

    /// Probability of a base matching the reference in a randomly mapped read
    double randomBaseMatchProb = 0.25;

    //------------------------------------------------------
    // continuous allele frequency caller options:

    // Assume a ploidy-based prior (0%, 50%, 100% or 0% 100% for haploid, diploid
    // If false, a continuous model is used
    bool is_ploidy_prior = true;

    /// The minimum allele frequency to call a heterozygous genotype when not assuming ploidy
    /// Must be in (0,0.5)
    double min_het_vf = 0.01;

    /// Expected allele observation quality used for the purpose of computing
    /// continuous frequency caller variant qualities
    const int continuousFrequencyCallerExpectedObservationQuality = 17;

    /// Maximum continuous frequency caller variant qscore
    const int continuousFrequencyCallerMaxQscore = 40;

    //------------------------------------------------------

    // vcfs can be input to specify candidate indels:
    std::vector<std::string> input_candidate_indel_vcf;

    // positions/indels in vcf must be written in output:
    std::vector<std::string> force_output_vcf;

    // Indicates that an upstream oligo is present on reads, which can be used to increase confidence for indels near the edge of the read
    unsigned upstream_oligo_size = 0;

    /// file specifying haploid and deleted regions as 1/0 in bed col 4
    std::string ploidy_region_vcf;

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

    bool useTier2Evidence = false;
    Tier2Options tier2;

    // indel error options
    std::vector<std::string> indelErrorModelFilenames;
    std::string thetaFilename;
    std::string indel_error_model_name = "logLinear";

    // Scalar multiple modifying the prob of observing an indel->reference error relative to reference->indel
    // This can be used to compensate for reference bias in the model.
    bool isIndelRefErrorFactor = false;
    double indelRefErrorFactor = 1.;

    /// When P(read | allele) is greater than or equal to this value, then the read counts as "supporting"
    /// the given allele
    ///
    /// WARNING: This value does not just change superficial AD count output in the VCF, but also impacts
    /// several count-based EVS metrics. An EVS retrain may be required when it is changed.
    PrettyFloat<double> readConfidentSupportThreshold = PrettyFloat<double>("0.51");

    // this option is only used by the error counting module
    //
    // per STREL-401, the error counting module excludes evidence close to the beginning or ending of each read,
    // so as to reduce systematic underestimation bias of indel error rates (or overestimation bias of SNV error rates)
    //
    unsigned minDistanceFromReadEdge = 0;

    /// Optional bedfile to specify which regions should be called in the genome
    std::string callRegionsBedFilename;

    /// If true, the original read alignment with soft-clipped edges is scored and chosen as the
    /// final alignment if it has the highest score.
    ///
    /// This should only be relevant when unrolling of soft-clipped read edges is used to improve indel calling,
    /// which is currently the case in Strelka.
    bool isRetainOptimalSoftClipping = false;
};



/// \brief Store sample-specific parameter values
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
struct GenotypePriorSet;


/// \brief Parameters deterministically derived from the input options
///
/// Note noncopyable required because of the unique ptrs.
struct starling_base_deriv_options : public blt_deriv_options, private boost::noncopyable
{
    typedef blt_deriv_options base_t;

    explicit
    starling_base_deriv_options(const starling_base_options& opt);

    ~starling_base_deriv_options();

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
    double randomBaseMatchLogProb;
    double tier2RandomBaseMatchLogProb;

    double correctMappingLogPrior;

    starling_align_limit sal;

    unsigned variant_window_first_stage;
    unsigned variant_window_last_stage;

    const min_count_binom_gte_cache countCache;

    // cache some frequently used log(rate) values
    const double logIndelRefErrorFactor;

private:
    std::unique_ptr<IndelErrorModel> _indelErrorModel;
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

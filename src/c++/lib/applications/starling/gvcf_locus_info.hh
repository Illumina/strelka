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

#include "gvcfAlleleInfo.hh"
#include "germlineVariantEmpiricalScoringFeatures.hh"
#include "ploidyUtil.hh"
#include "blt_common/position_snp_call_pprob_digt.hh"
#include "blt_util/align_path.hh"
#include "blt_util/math_util.hh"
#include "blt_util/PolymorphicObject.hh"
#include "starling_common/starling_indel_call_pprob_digt.hh"
#include "gvcf_options.hh"

#include <bitset>
#include <iosfwd>
#include <map>


namespace GERMLINE_VARIANT_VCF_FILTERS
{

enum index_t
{
    IndelConflict,
    SiteConflict,
    PloidyConflict,
    OffTarget,
    LowGQX,
    PhasingConflict,
    HighBaseFilt,
    HighDepth,
    HighSNVSB,
    HighSNVHPOL,
    HighRefRep,
    SIZE
};

inline
const char*
get_label(const unsigned idx)
{
    switch (idx)
    {
    case HighDepth:
        return "HighDepth";
    case LowGQX:
        return "LowGQX";
    case PhasingConflict:
        return "PhasingConflict";
    case HighSNVSB:
        return "HighSNVSB";
    case HighSNVHPOL:
        return "HighSNVHPOL";
    case HighBaseFilt:
        return "HighDPFRatio";
    case HighRefRep:
        return "HighREFREP";
    case IndelConflict:
        return "IndelConflict";
    case SiteConflict:
        return "SiteConflict";
    case PloidyConflict:
        return "PLOIDY_CONFLICT";
    case OffTarget:
        return "OffTarget";
    default:
        assert(false && "Unknown VCF filter value");
        return nullptr;
    }
}
}



struct GermlineFilterKeeper
{
    bool
    test(const GERMLINE_VARIANT_VCF_FILTERS::index_t i) const
    {
        return filters.test(i);
    }

    bool
    none() const
    {
        return filters.none();
    }

    bool
    any() const
    {
        return filters.any();
    }

    void
    set(const GERMLINE_VARIANT_VCF_FILTERS::index_t i)
    {
        filters.set(i);
    }

    void
    unset(const GERMLINE_VARIANT_VCF_FILTERS::index_t i)
    {
        filters.reset(i);
    }

    void
    write(std::ostream& os) const;

    bool
    operator==(const GermlineFilterKeeper& rhs) const
    {
        return (filters == rhs.filters);
    }

    // bit-wise or over each flag
    void
    merge(const GermlineFilterKeeper& filterKeeper)
    {
        filters |= filterKeeper.filters;
    }

    void
    clear()
    {
        filters.reset();
    }

private:
    std::bitset<GERMLINE_VARIANT_VCF_FILTERS::SIZE> filters;
};


/// property shared by all variants at a locus within one sample
struct LocusSampleInfo
{
    void
    clear()
    {
        ploidy.reset();
        gq = 0;
        filters.clear();
    }

    SamplePloidyState ploidy;
    int gq=0;
    GermlineFilterKeeper filters; ///< only for sample-specific filters
};


/// represents a locus in the sense of multiple alleles which interact in some way such that they would be represented in a single VCF record
struct LocusInfo : public PolymorphicObject
{
    explicit
    LocusInfo(
        const unsigned sampleCount,
        const pos_t initPos = 0)
      : pos(initPos),
       _sampleCount(sampleCount),
        _sampleInfo(sampleCount)
    {}

    void
    clear()
    {
        pos = 0;
        anyVariantAlleleQuality = 0;
        empiricalVariantScore = -1;
        filters.clear();
        for (auto& sample : _sampleInfo)
        {
            sample.clear();
        }
    }

    unsigned
    getSampleCount() const
    {
        return _sampleInfo.size();
    }

    LocusSampleInfo&
    getSample(const unsigned sampleIndex)
    {
        return _sampleInfo[sampleIndex];
    }

    const LocusSampleInfo&
    getSample(const unsigned sampleIndex) const
    {
        return _sampleInfo[sampleIndex];
    }

    /// zero-index position of the locus, alleles may not all start here:
    pos_t pos = 0;

    /// prob that there is a variant segregating at this sample (VCF calls this QUAL)
    int anyVariantAlleleQuality = 0;

    /// The empirically calibrated quality-score of the locus, if -1 no locus EVS is available
    int empiricalVariantScore = -1;

    /// All locus-ldevel filteres
    GermlineFilterKeeper filters;

private:
    unsigned _sampleCount;
    std::vector<LocusSampleInfo> _sampleInfo;
};


/// represents an indel call at the level of a full VCF record, containing possibly multiple alleles/SimpleGenotypes
struct GermlineIndelLocusInfo : public LocusInfo
{
    explicit
    GermlineIndelLocusInfo(
        const unsigned sampleCount,
        const pos_t initPos)
        : LocusInfo(sampleCount, initPos)
    {}

    virtual ~GermlineIndelLocusInfo() {}

    virtual bool isForcedOutput() const = 0;
    virtual bool is_indel() const = 0;
    // the EXCLUSIVE end of the variant (i.e. open)
    virtual pos_t end() const = 0;
};


/// specify that calling model is diploid
struct GermlineDiploidIndelLocusInfo : public GermlineIndelLocusInfo
{
    GermlineDiploidIndelLocusInfo(
        const gvcf_deriv_options& gvcfDerivedOptions,
        const unsigned sampleCount,
        const IndelKey& initIndelKey,
        const IndelData& initIndelData,
        const GermlineDiploidIndelSimpleGenotypeInfoCore& init_dindel,
        const starling_indel_report_info& initIndelReportInfo,
        const starling_indel_sample_report_info& initIndelSampleReportInfo)
        : GermlineIndelLocusInfo(sampleCount, initIndelKey.pos)
        , features(gvcfDerivedOptions.indelFeatureSet)
        , developmentFeatures(gvcfDerivedOptions.indelDevelopmentFeatureSet)

    {
        altAlleles.emplace_back(initIndelKey, initIndelData, initIndelReportInfo, initIndelSampleReportInfo, init_dindel);
    }

    bool isForcedOutput() const override
    {
        for (const auto& altAllele : altAlleles)
        {
            if (altAllele.isForcedOutput) return true;
        }
        return false;
    }

    bool is_indel() const override
    {
        return std::any_of(altAlleles.begin(), altAlleles.end(), [](const GermlineDiploidIndelAlleleInfo& x)
        {
            return x._dindel.isIndel();
        });
    }

    pos_t end() const override;

    void add_overlap(const reference_contig_segment& ref, GermlineDiploidIndelLocusInfo& overlap);

    ///TODO STREL-125 move this method to sample-level
    const char*
    get_gt() const
    {
        const auto& ploidy(getSample(0).ploidy);

        if (this->is_hetalt())
        {
            return "1/2";
        }
        else if (ploidy.isHaploid())
        {
            using namespace STAR_DIINDEL;

            switch (getFirstAltAllele().max_gt)
            {
            case NOINDEL:
                return "0";
            case HOM:
                return "1";
            default:
                assert(false && "Invalid indel genotype index");
                return "X";
            }
        }
        return STAR_DIINDEL::get_gt_label(getFirstAltAllele().max_gt);
    }

    ///TODO STREL-125 move this method to sample-level
    bool
    is_hetalt() const
    {
        return (_is_overlap);
    }

    ///TODO STREL-125 move this method to sample-level
    bool
    is_het() const
    {
        return (static_cast<int>(altAlleles.front().max_gt)>1);
    }

    // the site ploidy within the indel at offset x
    ///TODO STREL-125 move this method to sample-level
    unsigned
    getSitePloidy(const unsigned offset) const
    {
        const auto& ploidy(getSample(0).ploidy);

        auto& dindel(getFirstAltAllele()._dindel);
        if (ploidy.isNoploid()) return 0;

        if (!is_hetalt())
        {
            using namespace STAR_DIINDEL;
            switch (dindel.max_gt)
            {
            case HOM:
                return 0;
            case HET:
                return 1;
            case NOINDEL:
                return (ploidy.isHaploid() ? 1 : 2);
            }
            assert(0);
        }
        else
        {
            if (offset>=sitePloidy.size())
            {
                getPloidyError(offset);
            }
            return sitePloidy[offset];
        }
        return 2;
    }

    const GermlineDiploidIndelAlleleInfo& getFirstAltAllele() const
    {
        return altAlleles.front();
    }
    GermlineDiploidIndelAlleleInfo& getFirstAltAllele()
    {
        return altAlleles.front();
    }

    void
    computeEmpiricalScoringFeatures(
        const bool isRNA,
        const bool isUniformDepthExpected,
        const bool isComputeDevelopmentFeatures,
        const double chromDepth);

    void
    dump(std::ostream& os) const;

private:

    void
    getPloidyError(const unsigned offset) const;

public:

    /// represent site ploidy over the reference span of the overlapping indel set in the event of overlap:
    std::vector<unsigned> sitePloidy;

    // used to flag hetalt
    bool _is_overlap=false;

    std::vector<GermlineDiploidIndelAlleleInfo> altAlleles;

    /// production and development features used in the empirical scoring model:
    VariantScoringFeatureKeeper features;
    VariantScoringFeatureKeeper developmentFeatures;
};



/// represents an site call at the level of a full VCF record, containing possibly multiple alleles
struct GermlineSiteLocusInfo : public LocusInfo
{
    GermlineSiteLocusInfo(
        const unsigned sampleCount,
        const pos_t initPos,
        const char init_ref,
        const snp_pos_info& good_pi,
        const int used_allele_count_min_qscore,
        const bool initIsForcedOutput = false)
      : LocusInfo(sampleCount, initPos)
    {
        ref=(init_ref);
        isForcedOutput = initIsForcedOutput;
        good_pi.get_known_counts(fwd_counts,used_allele_count_min_qscore,true);
        good_pi.get_known_counts(rev_counts,used_allele_count_min_qscore,false);
        spanning_deletions = good_pi.n_spandel;
    }

    explicit
    GermlineSiteLocusInfo(
        const unsigned sampleCount)
        : LocusInfo(sampleCount)
    {}

    //GermlineSiteLocusInfo() = default;

    virtual bool is_snp() const = 0;
    virtual bool is_nonref() const = 0;

    unsigned
    alleleObservationCounts(const int base_id) const
    {
        return (fwd_counts[base_id]+rev_counts[base_id]);
    }

    unsigned
    alleleObservationCountsByStrand(
        const bool is_fwd_strand,
        const int base_id) const
    {
        return (is_fwd_strand ? fwd_counts[base_id] : rev_counts[base_id]);
    }

    void
    clear()
    {
        LocusInfo::clear();
        ref = 'N';
        n_used_calls = 0;
        n_unused_calls = 0;
        hpol = 0;
        spanning_deletions = 0;
        Unphasable = false;
        isForcedOutput = false;
    }

    char ref = 'N';
    unsigned n_used_calls = 0;
    unsigned n_unused_calls = 0;

    unsigned hpol = 0;

    unsigned spanning_deletions = 0;
    bool Unphasable = false;        // Set to true if the site should never be included in a phasing block
    bool isForcedOutput = false;

private:
    std::array<unsigned,N_BASE> fwd_counts;
    std::array<unsigned,N_BASE> rev_counts;
};


/// specify that calling model is diploid
struct GermlineDiploidSiteLocusInfo : public GermlineSiteLocusInfo
{
    GermlineDiploidSiteLocusInfo(
        const gvcf_deriv_options& gvcfDerivedOptions,
        const unsigned sampleCount,
        const pos_t init_pos,
        const char init_ref,
        const snp_pos_info& good_pi,
        const int used_allele_count_min_qscore,
        const bool is_forced_output = false)
        : GermlineSiteLocusInfo(sampleCount, init_pos, init_ref, good_pi, used_allele_count_min_qscore, is_forced_output),
          EVSFeatures(gvcfDerivedOptions.snvFeatureSet),
          EVSDevelopmentFeatures(gvcfDerivedOptions.snvDevelopmentFeatureSet)
    {}

    explicit
    GermlineDiploidSiteLocusInfo(
        const gvcf_deriv_options& gvcfDerivedOptions,
        const unsigned sampleCount)
        : GermlineSiteLocusInfo(sampleCount),
          EVSFeatures(gvcfDerivedOptions.snvFeatureSet),
          EVSDevelopmentFeatures(gvcfDerivedOptions.snvDevelopmentFeatureSet)
    {}

    bool is_snp() const override
    {
        return dgt.is_snp;
    }

    const char*
    get_gt() const
    {
        if       (allele.modified_gt != MODIFIED_SITE_GT::NONE)
        {
            return MODIFIED_SITE_GT::get_label(allele.modified_gt);
        }
        else if (is_print_unknowngt())
        {
            return ".";
        }
        else
        {
            const unsigned print_gt(allele.max_gt);
            return DIGT::get_vcf_gt(print_gt,dgt.ref_gt);
        }
    }

    void
    computeEmpiricalScoringFeatures(
        const bool isRNA,
        const bool isUniformDepthExpected,
        const bool isComputeDevelopmentFeatures,
        const double chromDepth);

    bool
    is_het() const
    {
        unsigned print_gt(allele.max_gt);
        return DIGT::is_het(print_gt);
    }

    bool
    is_hetalt() const
    {
        unsigned print_gt(allele.max_gt);
        const uint8_t a0(DIGT::get_allele(print_gt,0));
        const uint8_t a1(DIGT::get_allele(print_gt,1));
        return ((a0!=a1) && (dgt.ref_gt != a0) && (dgt.ref_gt != a1));
    }

    bool
    is_nonref() const override
    {
        return (allele.max_gt != dgt.ref_gt);
    }

    bool
    is_print_unknowngt() const
    {
        return (allele.is_unknown || (!allele.is_used_covered));
    }

    bool
    is_deletion() const
    {
        return ((!allele.is_unknown) && allele.is_used_covered && (!allele.is_zero_ploidy) && (is_nonref()));
    }

    bool
    is_qual() const
    {
        return ((!allele.is_unknown) && allele.is_used_covered && (!allele.is_zero_ploidy) && (is_nonref()));
    }

    void
    clearEVSFeatures()
    {
        EVSFeatures.clear();
        EVSDevelopmentFeatures.clear();
    }

    std::string phased_ref, phased_alt, phased_AD, phased_ADF, phased_ADR;
    diploid_genotype dgt;
    double hapscore = 0;
    double mapqRMS = 0;
    unsigned mapqZeroCount = 0;
    unsigned mapqCount = 0;

    //only meaningful for het calls
    double ReadPosRankSum = 0;  // Uses Mann-Whitney Rank Sum Test for the distance from the end of the read containing an alternate allele.
    double BaseQRankSum = 0;    // Uses Mann-Whitney Rank Sum Test for BQs (ref bases vs alternate alleles)
    double MQRankSum = 0;       // Uses Mann-Whitney Rank Sum Test for MQs (ref bases vs alternate alleles)
    double avgBaseQ = 0;
    double rawPos = 0;

    /// production and development features used in the empirical scoring model:
    VariantScoringFeatureKeeper EVSFeatures;
    VariantScoringFeatureKeeper EVSDevelopmentFeatures;

    GermlineDiploidSiteAlleleInfo allele;
};

std::ostream& operator<<(std::ostream& os,const GermlineDiploidSiteLocusInfo& si);


/// specify that variant is site/SNV and a continuous frequency calling model
struct GermlineContinuousSiteLocusInfo : public GermlineSiteLocusInfo
{
    GermlineContinuousSiteLocusInfo(
        const unsigned sampleCount,
        const pos_t init_pos,
        const char init_ref,
        const snp_pos_info& good_pi,
        const int used_allele_count_min_qscore,
        const double min_het_vf,
        const bool is_forced_output = false)
        : GermlineSiteLocusInfo(sampleCount,
                                init_pos,
                                init_ref,
                                good_pi,
                                used_allele_count_min_qscore,
                                is_forced_output),
          _min_het_vf(min_het_vf)
    {
    }

    bool is_snp() const override
    {
        return _is_snp;
    }

    bool is_nonref() const override
    {
        auto ref_id = base_to_id(ref);
        return altAlleles.end() !=
               std::find_if(altAlleles.begin(), altAlleles.end(),
                            [&](const GermlineContinuousSiteAlleleInfo& call)
        {
            return call._base != ref_id;
        });
    }

    const char* get_gt(const GermlineContinuousSiteAlleleInfo& call) const
    {
        if (call._base == base_to_id(ref))
            return "0/0";
        else if (call.variant_frequency() >= (1 -_min_het_vf))
            return "1/1";
        else if (call.variant_frequency() < _min_het_vf)
            return "0/0"; // STAR-66 - desired behavior
        else
            return "0/1";
    }

    bool _is_snp = false;
    std::vector<GermlineContinuousSiteAlleleInfo> altAlleles;
private:
    double _min_het_vf;
};


/// specify that variant is indel and a continuous frequency calling model
struct GermlineContinuousIndelLocusInfo : public GermlineIndelLocusInfo
{
    explicit
    GermlineContinuousIndelLocusInfo(
        const unsigned sampleCount,
        const pos_t init_pos)
        : GermlineIndelLocusInfo(sampleCount, init_pos)
    {
    }

    bool is_indel() const override
    {
        for (auto& call : altAlleles)
        {
            if (call._alleleDepth > 0)
                return true;
        }
        return false;
    }

    bool isForcedOutput() const override
    {
        for (const auto& altAllele : altAlleles)
        {
            if (altAllele.isForcedOutput) return true;
        }
        return false;
    }

    pos_t end() const override
    {
        pos_t result = 0;
        for (auto& x : altAlleles)
            result = std::max(result, x._indelKey.right_pos());
        return result;
    }

    const char* get_gt() const
    {
        if (is_het)
            return "0/1";
        else
            return "1/1";
    }

    std::vector<GermlineContinuousIndelAlleleInfo> altAlleles;
    bool is_het=false;
};

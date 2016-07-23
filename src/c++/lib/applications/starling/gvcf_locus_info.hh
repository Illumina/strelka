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

#include "germlineVariantEmpiricalScoringFeatures.hh"
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


/// this object contains information shared by any germline variant allele
///
/// variant here means SNV or indel
///
/// SimpleGenotype here means information pertaining to the core genotyping algorithm for SNVs or Indels (where for indels this genotype is
///    restricted to a single alt allele). This information is suplemnted with additional details requred for full VCF records in the "Call"
///    objects further below
///
/// model types where this is used include diploid/haploid and continuous single sample calling, but not contrastive (ie. tumor-normal) models
///
struct GermlineVariantSimpleGenotypeInfo : public PolymorphicObject
{
    GermlineVariantSimpleGenotypeInfo()
    {
        clear();
    }

    void
    set_filter(const GERMLINE_VARIANT_VCF_FILTERS::index_t i)
    {
        filters.set(i);
    }

    void
    unset_filter(const GERMLINE_VARIANT_VCF_FILTERS::index_t i)
    {
        filters.reset(i);
    }

    void
    write_filters(std::ostream& os) const;

    void
    clear()
    {
        gqx = 0;
        gq = 0;
        strand_bias = 0;
        empiricalVariantScore = -1;
        filters.reset();
    }

    int gqx=0;
    int gq=0;
    double strand_bias = 0;

    // The empirically calibrated quality-score of the site, if -1 no EVS is available
    int empiricalVariantScore = -1;

    std::bitset<GERMLINE_VARIANT_VCF_FILTERS::SIZE> filters;
};


std::ostream& operator<<(std::ostream& os,const GermlineVariantSimpleGenotypeInfo& shmod);


/// restrict to the case where variant is indel
struct GermlineIndelSimpleGenotypeInfo : public GermlineVariantSimpleGenotypeInfo
{
    GermlineIndelSimpleGenotypeInfo(
        const gvcf_deriv_options& gvcfDerivedOptions,
        const IndelKey& indelKey,
        const IndelData& indelData,
        const starling_indel_report_info indelReportInfo,
        const starling_indel_sample_report_info& indelSampleReportInfo)
        : _indelKey(indelKey)
        , _indelData(indelData)
        , _indelReportInfo(indelReportInfo)
        , _indelSampleReportInfo(indelSampleReportInfo)
        , features(gvcfDerivedOptions.indelFeatureSet)
        , developmentFeatures(gvcfDerivedOptions.indelDevelopmentFeatureSet)
    {}

    void
    clear()
    {
        GermlineVariantSimpleGenotypeInfo::clear();
        cigar.clear();
        features.clear();
        developmentFeatures.clear();
    }

    void set_hap_cigar(
        const unsigned lead=1,
        const unsigned trail=0);

    const IndelKey _indelKey;
    const IndelData _indelData;
    // TODO: make the indel overlapping code create a new call, then revert this to const
    starling_indel_report_info _indelReportInfo;
    const starling_indel_sample_report_info _indelSampleReportInfo;

    ALIGNPATH::path_t cigar;

    /// production and development features used in the empirical scoring model:
    VariantScoringFeatureKeeper features;
    VariantScoringFeatureKeeper developmentFeatures;
};

std::ostream& operator<<(std::ostream& os,const GermlineIndelSimpleGenotypeInfo& shi);


/// restrict to diploid calling models
struct GermlineDiploidIndelSimpleGenotypeInfo : public GermlineIndelSimpleGenotypeInfo
{
    GermlineDiploidIndelSimpleGenotypeInfo(
        const gvcf_deriv_options& gvcfDerivedOptions,
        const IndelKey& indelKey,
        const IndelData& indelData,
        const starling_indel_report_info& indelReportInfo,
        const starling_indel_sample_report_info& indelSampleReportInfo,
        const GermlineDiploidIndelSimpleGenotypeInfoCore& dindel)
        : GermlineIndelSimpleGenotypeInfo(gvcfDerivedOptions, indelKey, indelData, indelReportInfo, indelSampleReportInfo)
        , _dindel(dindel)
    {}

    void
    computeEmpiricalScoringFeatures(
        const bool isRNA,
        const bool isUniformDepthExpected,
        const bool isComputeDevelopmentFeatures,
        const double chromDepth,
        const bool isHetalt);

    /// TODO: this object is only here to trigger an auto copy ctor, seems like we can simplify now?
    /// TODO: Make indel_overlapper create new call objects, then revert this to const
    GermlineDiploidIndelSimpleGenotypeInfoCore _dindel;

    /// TODO: max_gt is here and in _dindel. Ugggghhh. Document the difference or get this down to one copy
    unsigned max_gt=0;
};

std::ostream& operator<<(std::ostream& os,const GermlineDiploidIndelSimpleGenotypeInfo& dic);



namespace MODIFIED_SITE_GT
{

enum index_t
{
    NONE,
    UNKNOWN,
    ZERO,
    ONE
};

inline
const char*
get_label(const unsigned idx)
{
    switch (static_cast<index_t>(idx))
    {
    case ZERO:
        return "0";
    case ONE:
        return "1";
    case UNKNOWN:
        return ".";
    default:
        assert(false && "Unknown site GT value");
        return nullptr;
    }
}
}


/// restrict to the case where variant is site/SNV and calling model is diploid
struct GermlineDiploidSiteSimpleGenotypeInfo : public GermlineVariantSimpleGenotypeInfo
{
    explicit
    GermlineDiploidSiteSimpleGenotypeInfo(
        const gvcf_deriv_options& gvcfDerivedOptions)
        : features(gvcfDerivedOptions.snvFeatureSet),
          developmentFeatures(gvcfDerivedOptions.snvDevelopmentFeatureSet)
    {
        clear();
    }

    void
    clear()
    {
        GermlineVariantSimpleGenotypeInfo::clear();
        is_unknown=true;
        is_covered=false;
        is_used_covered=false;
        is_zero_ploidy=false;
        is_phased_region=false;
        is_phasing_insufficient_depth=false;
        modified_gt=MODIFIED_SITE_GT::NONE;
        max_gt=0;
        clearEVSFeatures();
    }

    void
    clearEVSFeatures()
    {
        features.clear();
        developmentFeatures.clear();
    }

    bool
    is_gqx() const
    {
        return ((!is_unknown) && is_used_covered && (!is_zero_ploidy));
    }

    bool is_unknown;
    bool is_covered;
    bool is_used_covered;
    bool is_zero_ploidy;
    bool is_phased_region;
    bool is_phasing_insufficient_depth;

    MODIFIED_SITE_GT::index_t modified_gt;
    unsigned max_gt;

    /// production and development features used in the empirical scoring model:
    VariantScoringFeatureKeeper features;
    VariantScoringFeatureKeeper developmentFeatures;
};


/// restrict to the case where variant is site/SNV and calling model is continuous
struct GermlineContinuousSiteSimpleGenotypeInfo : public GermlineVariantSimpleGenotypeInfo
{
    GermlineContinuousSiteSimpleGenotypeInfo(unsigned totalDepth, unsigned alleleDepth, BASE_ID::index_t base)
        : _totalDepth(totalDepth)
        , _alleleDepth(alleleDepth)
        , _base(base)
    {
    }

    double variant_frequency() const
    {
        return safeFrac(_alleleDepth, _totalDepth);
    }

    unsigned _totalDepth;
    unsigned _alleleDepth;
    BASE_ID::index_t _base;
};


/// restrict to the case where variant is indel and calling model is continuous
struct GermlineContinuousIndelSimpleGenotypeInfo : public GermlineIndelSimpleGenotypeInfo
{
    GermlineContinuousIndelSimpleGenotypeInfo(
        const gvcf_deriv_options& gvcfDerivedOptions,
        unsigned totalDepth,
        unsigned alleleDepth,
        const IndelKey& indelKey,
        const IndelData& indelData,
        const starling_indel_report_info& indelReportInfo,
        const starling_indel_sample_report_info& indelSampleReportInfo)
        : GermlineIndelSimpleGenotypeInfo(gvcfDerivedOptions, indelKey, indelData, indelReportInfo, indelSampleReportInfo)
        , _totalDepth(totalDepth)
        , _alleleDepth(alleleDepth)
    {
        set_hap_cigar(0,0);
    }

    double variant_frequency() const
    {
        return safeFrac(_alleleDepth,_totalDepth);
    }

    unsigned _totalDepth;
    unsigned _alleleDepth;
};

std::ostream& operator<<(std::ostream& os,const GermlineDiploidSiteSimpleGenotypeInfo& smod);


/// represents an indel call at the level of a full VCF record, containing possibly multiple alleles/SimpleGenotypes
struct GermlineIndelCallInfo
{
    explicit GermlineIndelCallInfo(const pos_t init_pos)
        : pos(init_pos)
    {}

    virtual ~GermlineIndelCallInfo() {}

    virtual bool is_forced_output() const = 0;
    virtual bool is_indel() const = 0;
    virtual void set_filter(GERMLINE_VARIANT_VCF_FILTERS::index_t filter) = 0;
    // the EXCLUSIVE end of the variant (i.e. open)
    virtual pos_t end() const = 0;

    pos_t pos;
};


/// specify that calling model is diploid
struct GermlineDiploidIndelCallInfo : public GermlineIndelCallInfo
{
    GermlineDiploidIndelCallInfo(
        const gvcf_deriv_options& gvcfDerivedOptions,
        const IndelKey& initIndelKey,
        const IndelData& initIndelData,
        const GermlineDiploidIndelSimpleGenotypeInfoCore& init_dindel,
        const starling_indel_report_info& initIndelReportInfo,
        const starling_indel_sample_report_info& initIndelSampleReportInfo) : GermlineIndelCallInfo(initIndelKey.pos)
    {
        _calls.emplace_back(gvcfDerivedOptions, initIndelKey, initIndelData, initIndelReportInfo, initIndelSampleReportInfo, init_dindel);
    }

    bool is_forced_output() const override
    {
        return std::any_of(_calls.begin(), _calls.end(),
                           [](const GermlineDiploidIndelSimpleGenotypeInfo& x)
        {
            return x._dindel.is_forced_output;
        });
    }
    bool is_indel() const override
    {
        return std::any_of(_calls.begin(), _calls.end(), [](const GermlineDiploidIndelSimpleGenotypeInfo& x)
        {
            return x._dindel.isIndel();
        });
    }
    void set_filter(GERMLINE_VARIANT_VCF_FILTERS::index_t filter) override
    {
        for (auto& x : _calls) x.set_filter(filter);
    }

    pos_t end() const override;

    void add_overlap(const reference_contig_segment& ref, GermlineDiploidIndelCallInfo& overlap);

    const char*
    get_gt() const
    {
        if (this->is_hetalt())
        {
            return "1/2";
        }
        else if (first()._dindel.is_haploid())
        {
            using namespace STAR_DIINDEL;

            switch (first().max_gt)
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
        return STAR_DIINDEL::get_gt_label(first().max_gt);
    }

    bool
    is_hetalt() const
    {
        return (_is_overlap);
    }

    bool
    is_het() const
    {
        return (static_cast<int>(_calls.front().max_gt)>1);
    }

    // the site ploidy within the indel at offset x
    unsigned
    get_ploidy(const unsigned offset) const
    {
        auto& dindel(first()._dindel);
        if (dindel.is_noploid()) return 0;

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
                return (dindel.is_haploid() ? 1 : 2);
            }
            assert(0);
        }
        else
        {
            if (offset>=ploidy.size())
            {
                getPloidyError(offset);
            }
            return ploidy[offset];
        }
        return 2;
    }

    const GermlineDiploidIndelSimpleGenotypeInfo& first() const
    {
        return _calls.front();
    }
    GermlineDiploidIndelSimpleGenotypeInfo& first()
    {
        return _calls.front();
    }

    void
    dump(std::ostream& os) const;

private:

    void
    getPloidyError(const unsigned offset) const;

public:

    /// represent site ploidy over the reference span of the overlapping indel set in the event of overlap:
    std::vector<unsigned> ploidy;

    // used to flag hetalt
    bool _is_overlap=false;

    std::vector<GermlineDiploidIndelSimpleGenotypeInfo> _calls;
};



/// represents an site call at the level of a full VCF record, containing possibly multiple alleles
struct GermlineSiteCallInfo : public PolymorphicObject
{
    GermlineSiteCallInfo(const pos_t init_pos,
                         const char init_ref,
                         const snp_pos_info& good_pi,
                         const int used_allele_count_min_qscore,
                         const bool is_forced_output = false)
    {
        pos=(init_pos);
        ref=(init_ref);
        forcedOutput = is_forced_output;
        good_pi.get_known_counts(fwd_counts,used_allele_count_min_qscore,true);
        good_pi.get_known_counts(rev_counts,used_allele_count_min_qscore,false);
        spanning_deletions = good_pi.n_spandel;
    }

    GermlineSiteCallInfo() = default;

    virtual bool is_snp() const = 0;
    virtual void set_filter(GERMLINE_VARIANT_VCF_FILTERS::index_t filter) = 0;
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

    pos_t pos = 0;
    char ref = 'N';
    unsigned n_used_calls = 0;
    unsigned n_unused_calls = 0;

    unsigned hpol = 0;

    unsigned spanning_deletions;
    bool Unphasable = false;        // Set to true if the site should never be included in a phasing block
    bool forcedOutput = false;

private:
    std::array<unsigned,N_BASE> fwd_counts;
    std::array<unsigned,N_BASE> rev_counts;
};


/// specify that calling model is diploid
struct GermlineDiploidSiteCallInfo : public GermlineSiteCallInfo
{
    GermlineDiploidSiteCallInfo(
        const gvcf_deriv_options& gvcfDerivedOptions,
        const pos_t init_pos,
        const char init_ref,
        const snp_pos_info& good_pi,
        const int used_allele_count_min_qscore,
        const bool is_forced_output = false)
        : GermlineSiteCallInfo(init_pos, init_ref, good_pi, used_allele_count_min_qscore, is_forced_output),
          smod(gvcfDerivedOptions)
    {}

    explicit
    GermlineDiploidSiteCallInfo(
        const gvcf_deriv_options& gvcfDerivedOptions)
        : smod(gvcfDerivedOptions)
    {}

    bool is_snp() const override
    {
        return dgt.is_snp;
    }

    void set_filter (GERMLINE_VARIANT_VCF_FILTERS::index_t filter) override
    {
        smod.set_filter(filter);
    }

    const char*
    get_gt() const
    {
        if       (smod.modified_gt != MODIFIED_SITE_GT::NONE)
        {
            return MODIFIED_SITE_GT::get_label(smod.modified_gt);
        }
        else if (is_print_unknowngt())
        {
            return ".";
        }
        else
        {
            const unsigned print_gt(smod.max_gt);
            return DIGT::get_vcf_gt(print_gt,dgt.ref_gt);
        }
    }

    void
    computeEmpiricalScoringFeatures(
        const bool isRNA,
        const bool isUniformDepthExpected,
        const bool isComputeDevelopmentFeatures,
        const double chromDepth,
        GermlineDiploidSiteSimpleGenotypeInfo& smod2) const;

    bool
    is_het() const
    {
        unsigned print_gt(smod.max_gt);
        return DIGT::is_het(print_gt);
    }

    bool
    is_hetalt() const
    {
        unsigned print_gt(smod.max_gt);
        const uint8_t a0(DIGT::get_allele(print_gt,0));
        const uint8_t a1(DIGT::get_allele(print_gt,1));
        return ((a0!=a1) && (dgt.ref_gt != a0) && (dgt.ref_gt != a1));
    }

    bool
    is_nonref() const override
    {
        return (smod.max_gt != dgt.ref_gt);
    }

    bool
    is_print_unknowngt() const
    {
        return (smod.is_unknown || (!smod.is_used_covered));
    }

    bool
    is_deletion() const
    {
        return ((!smod.is_unknown) && smod.is_used_covered && (!smod.is_zero_ploidy) && (is_nonref()));
    }

    bool
    is_qual() const
    {
        return ((!smod.is_unknown) && smod.is_used_covered && (!smod.is_zero_ploidy) && (is_nonref()));
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

    GermlineDiploidSiteSimpleGenotypeInfo smod;
};

std::ostream& operator<<(std::ostream& os,const GermlineDiploidSiteCallInfo& si);


/// specify that variant is site/SNV and a continuous frequency calling model
struct GermlineContinuousSiteCallInfo : public GermlineSiteCallInfo
{
    GermlineContinuousSiteCallInfo(
        const pos_t init_pos,
        const char init_ref,
        const snp_pos_info& good_pi,
        const int used_allele_count_min_qscore,
        const double min_het_vf,
        const bool is_forced_output = false) : GermlineSiteCallInfo(init_pos,
                                                                        init_ref,
                                                                        good_pi,
                                                                        used_allele_count_min_qscore,
                                                                        is_forced_output)
        , _min_het_vf(min_het_vf)
    {
    }

    bool is_snp() const override
    {
        return _is_snp;
    }
    void set_filter (GERMLINE_VARIANT_VCF_FILTERS::index_t filter) override
    {
        for (auto& call : calls) call.set_filter(filter);
    }
    bool is_nonref() const override
    {
        auto ref_id = base_to_id(ref);
        return calls.end() !=
               std::find_if(calls.begin(), calls.end(),
                            [&](const GermlineContinuousSiteSimpleGenotypeInfo& call)
        {
            return call._base != ref_id;
        });
    }

    const char* get_gt(const GermlineContinuousSiteSimpleGenotypeInfo& call) const
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
    std::vector<GermlineContinuousSiteSimpleGenotypeInfo> calls;
private:
    double _min_het_vf;
};


/// specify that variant is indel and a continuous frequency calling model
struct GermlineContinuousIndelCallInfo : public GermlineIndelCallInfo
{
    explicit GermlineContinuousIndelCallInfo(const pos_t init_pos)
        : GermlineIndelCallInfo(init_pos)
    {
    }

    void set_filter(GERMLINE_VARIANT_VCF_FILTERS::index_t filter) override
    {
        for (auto& call : calls) call.set_filter(filter);
    }


    bool is_indel() const override
    {
        for (auto& call : calls)
        {
            if (call._alleleDepth > 0)
                return true;
        }
        return false;
    }

    bool is_forced_output() const override
    {
        for (auto& call : calls)
            if (call._indelData.is_forced_output)
                return true;
        return false;
    }

    pos_t end() const override
    {
        pos_t result = 0;
        for (auto& x : calls)
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

    std::vector<GermlineContinuousIndelSimpleGenotypeInfo> calls;
    bool is_het=false;
};




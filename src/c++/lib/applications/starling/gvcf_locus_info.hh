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



struct LocusFilterKeeper
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
    operator==(const LocusFilterKeeper& rhs) const
    {
        return (filters == rhs.filters);
    }

    // bit-wise or over each flag
    void
    merge(const LocusFilterKeeper& filterKeeper)
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



/// represents a locus in the sense of multiple alleles which interact in some way such that they would be represented in a single VCF record
struct LocusInfo : public PolymorphicObject
{
    void
    clear()
    {
        filters.clear();
    }

    LocusFilterKeeper filters;
};


/// represents an indel call at the level of a full VCF record, containing possibly multiple alleles/SimpleGenotypes
struct GermlineIndelLocusInfo : public LocusInfo
{
    explicit GermlineIndelLocusInfo(const pos_t init_pos)
        : pos(init_pos)
    {}

    virtual ~GermlineIndelLocusInfo() {}

    virtual bool is_forced_output() const = 0;
    virtual bool is_indel() const = 0;
    // the EXCLUSIVE end of the variant (i.e. open)
    virtual pos_t end() const = 0;

    pos_t pos;
};


/// specify that calling model is diploid
struct GermlineDiploidIndelLocusInfo : public GermlineIndelLocusInfo
{
    GermlineDiploidIndelLocusInfo(
        const gvcf_deriv_options& gvcfDerivedOptions,
        const IndelKey& initIndelKey,
        const IndelData& initIndelData,
        const GermlineDiploidIndelSimpleGenotypeInfoCore& init_dindel,
        const starling_indel_report_info& initIndelReportInfo,
        const starling_indel_sample_report_info& initIndelSampleReportInfo) : GermlineIndelLocusInfo(initIndelKey.pos)
    {
        altAlleles.emplace_back(gvcfDerivedOptions, initIndelKey, initIndelData, initIndelReportInfo, initIndelSampleReportInfo, init_dindel);
    }

    bool is_forced_output() const override
    {
        return std::any_of(altAlleles.begin(), altAlleles.end(),
                           [](const GermlineDiploidIndelAlleleInfo& x)
        {
            return x._dindel.is_forced_output;
        });
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
        return (static_cast<int>(altAlleles.front().max_gt)>1);
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

    const GermlineDiploidIndelAlleleInfo& first() const
    {
        return altAlleles.front();
    }
    GermlineDiploidIndelAlleleInfo& first()
    {
        return altAlleles.front();
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

    std::vector<GermlineDiploidIndelAlleleInfo> altAlleles;
};



/// represents an site call at the level of a full VCF record, containing possibly multiple alleles
struct GermlineSiteLocusInfo : public LocusInfo
{
    GermlineSiteLocusInfo(
        const pos_t init_pos,
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

    GermlineSiteLocusInfo() = default;

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
        pos = 0;
        ref = 'N';
        n_used_calls = 0;
        n_unused_calls = 0;
        hpol = 0;
        spanning_deletions = 0;
        Unphasable = false;
        forcedOutput = false;
    }

    pos_t pos = 0;
    char ref = 'N';
    unsigned n_used_calls = 0;
    unsigned n_unused_calls = 0;

    unsigned hpol = 0;

    unsigned spanning_deletions = 0;
    bool Unphasable = false;        // Set to true if the site should never be included in a phasing block
    bool forcedOutput = false;

private:
    std::array<unsigned,N_BASE> fwd_counts;
    std::array<unsigned,N_BASE> rev_counts;
};


/// specify that calling model is diploid
struct GermlineDiploidSiteLocusInfo : public GermlineSiteLocusInfo
{
    GermlineDiploidSiteLocusInfo(
        const gvcf_deriv_options& gvcfDerivedOptions,
        const pos_t init_pos,
        const char init_ref,
        const snp_pos_info& good_pi,
        const int used_allele_count_min_qscore,
        const bool is_forced_output = false)
        : GermlineSiteLocusInfo(init_pos, init_ref, good_pi, used_allele_count_min_qscore, is_forced_output),
          EVSFeatures(gvcfDerivedOptions.snvFeatureSet),
          EVSDevelopmentFeatures(gvcfDerivedOptions.snvDevelopmentFeatureSet)
    {}

    explicit
    GermlineDiploidSiteLocusInfo(
        const gvcf_deriv_options& gvcfDerivedOptions)
        : EVSFeatures(gvcfDerivedOptions.snvFeatureSet),
          EVSDevelopmentFeatures(gvcfDerivedOptions.snvDevelopmentFeatureSet)
    {}

    bool is_snp() const override
    {
        return dgt.is_snp;
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
        const double chromDepth);

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

    GermlineDiploidSiteAlleleInfo smod;
};

std::ostream& operator<<(std::ostream& os,const GermlineDiploidSiteLocusInfo& si);


/// specify that variant is site/SNV and a continuous frequency calling model
struct GermlineContinuousSiteLocusInfo : public GermlineSiteLocusInfo
{
    GermlineContinuousSiteLocusInfo(
        const pos_t init_pos,
        const char init_ref,
        const snp_pos_info& good_pi,
        const int used_allele_count_min_qscore,
        const double min_het_vf,
        const bool is_forced_output = false) : GermlineSiteLocusInfo(init_pos,
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
    explicit GermlineContinuousIndelLocusInfo(const pos_t init_pos)
        : GermlineIndelLocusInfo(init_pos)
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

    bool is_forced_output() const override
    {
        for (auto& call : altAlleles)
            if (call._indelData.is_forced_output)
                return true;
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

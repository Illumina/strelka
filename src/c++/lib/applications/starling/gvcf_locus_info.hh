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
#include "gvcf_options.hh"
#include "germlineVariantEmpiricalScoringFeatures.hh"
#include "ploidyUtil.hh"
#include "blt_common/position_snp_call_pprob_digt.hh"
#include "blt_util/align_path.hh"
#include "blt_util/math_util.hh"
#include "blt_util/PolymorphicObject.hh"
#include "htsapi/vcf_util.hh"
#include "starling_common/starling_indel_call_pprob_digt.hh"
#include "starling_common/LocusSupportingReadStats.hh"

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

    // bit-wise and over each flag
    void
    unionMerge(const GermlineFilterKeeper& filterKeeper)
    {
        filters &= filterKeeper.filters;
    }

    void
    clear()
    {
        filters.reset();
    }

private:
    std::bitset<GERMLINE_VARIANT_VCF_FILTERS::SIZE> filters;
};


///
///
/// store data indexed by genotype, as ordered by the VCF standard
///
struct GenotypeLikelihoods
{
    typedef std::vector<unsigned> data_t;
    typedef data_t::const_iterator const_iterator;

    const_iterator
    begin() const
    {
        return _genotypePhredLoghood.begin();
    }

    const_iterator
    end() const
    {
        return _genotypePhredLoghood.end();
    }

    void
    clear()
    {
        _ploidy=0;
        _genotypePhredLoghood.clear();
    }

    void
    setPloidy(const int ploidy)
    {
        _ploidy=ploidy;
    }

    std::vector<unsigned>&
    getGenotypeLikelihood()
    {
        return _genotypePhredLoghood;
    }

    const std::vector<unsigned>&
    getGenotypeLikelihood() const
    {
        return _genotypePhredLoghood;
    }

    unsigned&
    getGenotypeLikelihood(
        const uint8_t allele0Index)
    {
        const unsigned index(getGenotypeIndex(allele0Index));
        if (index>=_genotypePhredLoghood.size())
        {
            _genotypePhredLoghood.resize(index+1,0);
        }
        return _genotypePhredLoghood[index];
    }

    unsigned
    getGenotypeLikelihood(
        const uint8_t allele0Index) const
    {
        const unsigned index(getGenotypeIndex(allele0Index));
        assert(index<_genotypePhredLoghood.size());
        return _genotypePhredLoghood[index];
    }

    unsigned&
    getGenotypeLikelihood(
        const uint8_t allele0Index,
        const uint8_t allele1Index)
    {
        const unsigned index(getGenotypeIndex(allele0Index, allele1Index));
        if (index>=_genotypePhredLoghood.size())
        {
            const unsigned index2(getGenotypeIndex(allele1Index, allele1Index));
            _genotypePhredLoghood.resize(index2+1,0);
        }
        return _genotypePhredLoghood[index];
    }

    unsigned
    getGenotypeLikelihood(
        const uint8_t allele0Index,
        const uint8_t allele1Index) const
    {
        const unsigned index(getGenotypeIndex(allele0Index, allele1Index));
        assert(index<_genotypePhredLoghood.size());
        return _genotypePhredLoghood[index];
    }

private:
    unsigned
    getGenotypeIndex(
        const uint8_t allele0Index) const
    {
        assert(_ploidy==1);
        return VcfGenotypeUtil::getGenotypeIndex(allele0Index);
    }

    unsigned
    getGenotypeIndex(
        const uint8_t allele0Index,
        const uint8_t allele1Index) const
    {
        assert(_ploidy==2);
        return VcfGenotypeUtil::getGenotypeIndex(allele0Index,allele1Index);
    }

    int _ploidy = 0;

    /// likelihoods for all possible genotypes, defined as a function of ploidy and altAllele count
    std::vector<unsigned> _genotypePhredLoghood;
};



/// property shared by all variants at a locus within one sample
struct LocusSampleInfo
{
    void
    clear()
    {
        genotypeQualityPolymorphic = 0;
        maxGenotypeIndexPolymorphic = 0;
        genotypeQuality = 0;
        maxGenotypeIndex = 0;
        gqx = 0;
        empiricalVariantScore = -1;
        genotypePhredLoghood.clear();
        filters.clear();
        supportCounts.clear();
        _ploidy.reset();
        _isPloidyConflict = false;
    }

    const SamplePloidyState&
    getPloidy() const
    {
        return _ploidy;
    }

    void
    setPloidy(const int ploidy)
    {
        _ploidy.setPloidy(ploidy);
        genotypePhredLoghood.setPloidy(ploidy);
    }

    bool
    isPloidyConflict() const
    {
        return _isPloidyConflict;
    }

    void
    setPloidyConflict()
    {
        _isPloidyConflict=true;
    }

    unsigned max_gt() const
    {
        return maxGenotypeIndexPolymorphic;
    }

    /// non-forced sample printing criteria
    /// (for indels at least)
    bool
    isVariant() const
    {
        return (max_gt() != 0);
    }

    //--------------------------------------------
    // data

    /// VCF GQ
    int genotypeQualityPolymorphic=0;

    /// VCF GT
    unsigned maxGenotypeIndexPolymorphic=0;

    int genotypeQuality=0;
    unsigned maxGenotypeIndex=0;

    /// VCF GQX
    int gqx=0;

    /// likelihoods for all possible genotypes, defined as a function of ploidy and altAllele count
    /// used to write VCF PL vals
    GenotypeLikelihoods genotypePhredLoghood;

    /// The empirically calibrated quality-score of the sample genotype, if -1 no locus EVS is available
    int empiricalVariantScore = -1;

    /// only for sample-specific filters
    GermlineFilterKeeper filters;

    /// VCF AD/ADF/ADR counts
    LocusSupportingReadStats supportCounts;

private:
    SamplePloidyState _ploidy;
    bool _isPloidyConflict = false;
};

std::ostream& operator<<(std::ostream& os,const LocusSampleInfo& lsi);


/// represents a locus in the sense of multiple alleles which interact in some way such that they would be represented in a single VCF record
struct LocusInfo : public PolymorphicObject
{
    explicit
    LocusInfo(
        const unsigned sampleCount,
        const pos_t initPos = 0)
      : pos(initPos),
        _sampleInfo(sampleCount)
    {}

    void
    clear()
    {
        pos = 0;
        anyVariantAlleleQuality = 0;
        filters.clear();
        for (auto& sample : _sampleInfo)
        {
            sample.clear();
        }
        _altAlleleCount=0;
    }

    unsigned
    getAltAlleleCount() const
    {
        return _altAlleleCount;
    }

    void
    incrementAltAlleleCount()
    {
        _altAlleleCount++;
    }

    unsigned
    getSampleCount() const
    {
        return _sampleInfo.size();
    }

    LocusSampleInfo&
    getSample(const unsigned sampleIndex)
    {
        //assert(sampleIndex<_sampleInfo.size());
        return _sampleInfo[sampleIndex];
    }

    const LocusSampleInfo&
    getSample(const unsigned sampleIndex) const
    {
        //assert(sampleIndex<_sampleInfo.size());
        return _sampleInfo[sampleIndex];
    }

    /// non-forced locus printing criteria
    // (for diploid indels at least)
    bool
    isVariantLocus() const
    {
        for (const auto& sample : _sampleInfo)
        {
            if (sample.isVariant()) return true;
        }
        return false;
    }

    /// zero-index position of the locus, alleles may not all start here:
    pos_t pos = 0;

    /// prob that any of the ALTs exist in any of the samples (VCF calls this QUAL)
    int anyVariantAlleleQuality = 0;

    /// All locus-level filters
    GermlineFilterKeeper filters;

private:
    std::vector<LocusSampleInfo> _sampleInfo;
    unsigned _altAlleleCount = 0;
};

std::ostream& operator<<(std::ostream& os,const LocusInfo& li);


struct GermlineIndelSampleInfo
{
    double alleleFrequency() const
    {
        return safeFrac(reportInfo.n_confident_indel_reads, reportInfo.total_confident_reads());
    }

    AlleleSampleReportInfo reportInfo;

    /// the expected ploidy of sites spanning the indel locus assuming the ML GT is true
    std::vector<uint8_t> sitePloidy;
};


/// represents an indel call at the level of a full VCF record, containing possibly multiple alleles/SimpleGenotypes
struct GermlineIndelLocusInfo : public LocusInfo
{
    explicit
    GermlineIndelLocusInfo(
        const unsigned sampleCount)
        : LocusInfo(sampleCount),
          _indelSampleInfo(sampleCount)
    {}

    virtual ~GermlineIndelLocusInfo() {}

    const known_pos_range2&
    range() const
    {
        assert(getAltAlleleCount() > 0);
        return _range;
    }

    pos_t
    end() const
    {
        assert(getAltAlleleCount() > 0);
        return _range.end_pos();
    }

    void
    setIndelSampleInfo(
        const unsigned sampleIndex,
        const GermlineIndelSampleInfo& indelSampleInfo)
    {
        // ensure that no alleles are added once we start adding samples...
        assert(getAltAlleleCount()>0);
        _isLockAlleles = true;
        assert(indelSampleInfo.sitePloidy.size() == _range.size());
        assert(sampleIndex < _indelSampleInfo.size());
        _indelSampleInfo[sampleIndex] = indelSampleInfo;
    }

    const GermlineIndelSampleInfo&
    getIndelSample(const unsigned sampleIndex) const
    {
        return _indelSampleInfo[sampleIndex];
    }

    void
    addAltIndelAllele(
        const IndelKey& indelKey,
        const IndelData& indelData)
    {
        assert(not _isLockAlleles);

        _indelAlleleInfo.emplace_back(indelKey, indelData);
        incrementAltAlleleCount();
        if (_indelAlleleInfo.size() == 1)
        {
            _range.set_range(indelKey.pos, indelKey.right_pos());
            pos=_range.begin_pos();
        }
        else
        {
            _range.merge_range(known_pos_range2(indelKey.pos, indelKey.right_pos()));
            if (pos > _range.begin_pos()) pos = _range.begin_pos();
        }


    }

    const std::vector<GermlineIndelAlleleInfo>&
    getIndelAlleles() const
    {
        return _indelAlleleInfo;
    }

    /// get the ploidy of a site spanned by the maxGT indels at this locus in specified sample:
    unsigned
    getSitePloidy(
        const unsigned sampleIndex,
        const unsigned offset) const
    {
        if (offset>=_range.size())
        {
            getOffsetError(offset);
        }
        const GermlineIndelSampleInfo& indelSampleInfo(getIndelSample(sampleIndex));
        assert(_range.size() == indelSampleInfo.sitePloidy.size());
        return indelSampleInfo.sitePloidy[offset];
    }

    bool
    isAnyForcedOutputAtLocus() const
    {
        for (const auto& allele : _indelAlleleInfo)
        {
            if (allele.isForcedOutput) return true;
        }
        return false;
    }

    bool
    isAnyBreakpointAlleles() const
    {
        for (const auto& allele : _indelAlleleInfo)
        {
            if (allele.indelKey.is_breakpoint()) return true;
        }
        return false;
    }

    /// run all internal consistency checks
    void
    assertValidity() const
    {
        const unsigned sampleCount(getSampleCount());
        assert (sampleCount == _indelSampleInfo.size());
        for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
        {
            const GermlineIndelSampleInfo& indelSampleInfo(getIndelSample(sampleIndex));
            assert(_range.size() == indelSampleInfo.sitePloidy.size());
        }
    }

private:

    void
    getOffsetError(const unsigned offset) const;

    std::vector<GermlineIndelAlleleInfo> _indelAlleleInfo;
    std::vector<GermlineIndelSampleInfo> _indelSampleInfo;

    /// the refernece range of all indel alleles at this locus:
    known_pos_range2 _range;

    /// to sanity check input, locus must be specified by adding all alleles, and then adding all sample information, this bool enforces the allele->sample ordering
    bool _isLockAlleles = false;
};

std::ostream& operator<<(std::ostream& os,const GermlineIndelLocusInfo& ii);


/// specify that calling model is diploid
struct GermlineDiploidIndelLocusInfo : public GermlineIndelLocusInfo
{
    GermlineDiploidIndelLocusInfo(
        const gvcf_deriv_options& gvcfDerivedOptions,
        const unsigned sampleCount)
        : GermlineIndelLocusInfo(sampleCount)
        , evsFeatures(gvcfDerivedOptions.indelFeatureSet)
        , evsDevelopmentFeatures(gvcfDerivedOptions.indelDevelopmentFeatureSet)
    {}

    /// \param allSampleChromDepth expected depth summed over all samples
    static
    void
    computeEmpiricalScoringFeatures(
        const GermlineDiploidIndelLocusInfo& locus,
        const unsigned sampleIndex,
        const bool isRNA,
        const bool isUniformDepthExpected,
        const bool isComputeDevelopmentFeatures,
        const double allSampleChromDepth,
        VariantScoringFeatureKeeper& features,
        VariantScoringFeatureKeeper& developmentFeatures);

    void
    dump(std::ostream& os) const;

    /// production and development features used in the empirical scoring model:
    VariantScoringFeatureKeeper evsFeatures;
    VariantScoringFeatureKeeper evsDevelopmentFeatures;
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

std::ostream& operator<<(std::ostream& os,const GermlineSiteLocusInfo& si);


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
          evsFeatures(gvcfDerivedOptions.snvFeatureSet),
          evsDevelopmentFeatures(gvcfDerivedOptions.snvDevelopmentFeatureSet)
    {}

    explicit
    GermlineDiploidSiteLocusInfo(
        const gvcf_deriv_options& gvcfDerivedOptions,
        const unsigned sampleCount)
        : GermlineSiteLocusInfo(sampleCount),
          evsFeatures(gvcfDerivedOptions.snvFeatureSet),
          evsDevelopmentFeatures(gvcfDerivedOptions.snvDevelopmentFeatureSet)
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

    /// \param allSampleChromDepth expected depth summed over all samples
    static
    void
    computeEmpiricalScoringFeatures(
        const GermlineDiploidSiteLocusInfo& locus,
        const unsigned sampleIndex,
        const bool isRNA,
        const bool isUniformDepthExpected,
        const bool isComputeDevelopmentFeatures,
        const double allSampleChromDepth,
        VariantScoringFeatureKeeper& features,
        VariantScoringFeatureKeeper& developmentFeatures);

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
        evsFeatures.clear();
        evsDevelopmentFeatures.clear();
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
    VariantScoringFeatureKeeper evsFeatures;
    VariantScoringFeatureKeeper evsDevelopmentFeatures;

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
        const unsigned sampleCount)
        : GermlineIndelLocusInfo(sampleCount)
    {}
};

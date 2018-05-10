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

/// \file
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
    LowGQX,
    PhasingConflict,
    HighBaseFilt,
    HighDepth,
    HighSNVSB,
    HighSNVHPOL,
    HighRefRep,
    LowDepth,
    NoPassedVariantGTs,
    NotGenotyped,
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
        return "PloidyConflict";
    case LowDepth:
        return "LowDepth";
    case NoPassedVariantGTs:
        return "NoPassedVariantGTs";
    case NotGenotyped:
        return "NotGenotyped";
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
    intersectWith(const GermlineFilterKeeper& filterKeeper)
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
        maxGenotypeIndexPolymorphic.clear();
        genotypeQuality = 0;
        maxGenotypeIndex.clear();
        gqx = 0;
        empiricalVariantScore = -1;
        genotypePhredLoghood.clear();
        filters.clear();
        supportCounts.clear();
        _activeRegionId = -1;
        phaseSetId = -1;
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

    VcfGenotype&
    max_gt()
    {
        return maxGenotypeIndexPolymorphic;
    }

    const VcfGenotype&
    max_gt() const
    {
        return maxGenotypeIndexPolymorphic;
    }

    /// Return true if the MAP genotype includes a non-reference allele.
    ///
    /// This is currently the non-forced sample printing criteria (for indels at least)
    bool
    isVariant() const
    {
        return (max_gt().isVariant());
    }

    /// sets gqx given that other relevant values have been setup:
    void
    setGqx()
    {
        // maxGenotypeIndex != maxGenotypeIndexPolymorphic indicates we're in a boundary zone
        // between variant and hom-ref call
        if (not (maxGenotypeIndex == maxGenotypeIndexPolymorphic))
        {
            gqx = 0;
        }
        else
        {
            gqx = std::min(genotypeQuality, genotypeQualityPolymorphic);
        }
    }

    void
    setActiveRegionId(ActiveRegionId activeRegionId)
    {
        _activeRegionId = activeRegionId;
    }

    ActiveRegionId
    getActiveRegionId() const
    {
        return _activeRegionId;
    }

    //--------------------------------------------
    // data

    /// VCF GQ
    int genotypeQualityPolymorphic=0;

    /// VCF GT
    VcfGenotype maxGenotypeIndexPolymorphic;

    int genotypeQuality=0;
    VcfGenotype maxGenotypeIndex;

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


    /// if non-negative, use this to provide a PS tag to accompany a phased max_gt() value.
    pos_t phaseSetId = -1;

private:
    SamplePloidyState _ploidy;
    bool _isPloidyConflict = false;

    /// Active region ID (-1 if this locus is not in AR)
    ActiveRegionId _activeRegionId = -1;

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

    /// locus is variant if any sample has an ALT in its most likely GT
    ///
    /// previous designs used QUAL>0, consider using both?
    bool
    isVariantLocus() const
    {
        for (const auto& sample : _sampleInfo)
        {
            if (sample.isVariant()) return true;
        }
        return false;
    }

    /// \return Minimum active region ID among all samples.
    /// The return value is -1 if this locus is not in AR in any of the samples.
    /// The maximum valid active region ID is the current position.
    ActiveRegionId getActiveRegionId() const
    {
        ActiveRegionId activeRegionId(pos+1);

        for (const auto& sample : _sampleInfo)
        {
            auto sampleActiveRegionId(sample.getActiveRegionId());
            if ((sampleActiveRegionId >= 0) && (sampleActiveRegionId < activeRegionId))
                activeRegionId = sampleActiveRegionId;
        }

        if (activeRegionId > pos)
            return -1;  // this locus is not in AR in any of the samples

        return activeRegionId;
    }

protected:
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
        _isLockAlleles = false;
    }

public:
    // ---------------------------------------------------------
    /// zero-index position of the locus, alleles may not all start here:
    pos_t pos = 0;

    /// prob that any of the ALTs exist in any of the samples (VCF calls this QUAL)
    int anyVariantAlleleQuality = 0;

    /// All locus-level filters
    GermlineFilterKeeper filters;

private:
    std::vector<LocusSampleInfo> _sampleInfo;
    unsigned _altAlleleCount = 0;

protected:
    /// to sanity check input, locus must be specified by adding all alleles, and then adding all sample information, this bool enforces the allele->sample ordering
    bool _isLockAlleles = false;

};

std::ostream& operator<<(std::ostream& os,const LocusInfo& li);


struct GermlineIndelSampleInfo
{
    /// get the count of reads crossing indel or reference alleles at any teir or mapq value, including MAPQ0 reads
    ///
    /// TODO get a more exact definition of requirements for read to be counted here - how are mulitple indel alleles handled, by how many bases must the read cross an indel breakend? etc...
    unsigned
    getTotalReadDepth() const
    {
        return mapqTracker.count;
    }

    /// the depth of the pileup at the position preceding the locus
    unsigned tier1Depth = 0;

    /// mapq stats at the position preceding the locus
    MapqTracker mapqTracker;


    /// TODO STREL-125 deprecated - replaced by supportCounts in LocusSampleInfo
    double alleleFrequency() const
    {
        return safeFrac(legacyReportInfo.n_confident_indel_reads, legacyReportInfo.total_confident_reads());
    }

    /// TODO STREL-125 deprecated - replaced by supportCounts in LocusSampleInfo
    AlleleSampleReportInfo legacyReportInfo;
};


/// represents an indel call at the level of a full VCF record, containing possibly multiple alleles/SimpleGenotypes
struct GermlineIndelLocusInfo : public LocusInfo
{
    explicit
    GermlineIndelLocusInfo(
        const unsigned sampleCount)
        : LocusInfo(sampleCount),
          _indelSampleInfo(sampleCount), _commonPrefixLength(0), _doNotGenotype(false)
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

    /// Fix locus for the case that all alleles share a common prefix
    ///
    /// This was added to address STREL-275. Ideally we would design the genotyping method to never input
    /// shared prefix loci in the first place, This is a temporary patch to prevent invalid output for now.
    void
    setCommonPrefix(const unsigned commonPrefixLength)
    {
        _commonPrefixLength = commonPrefixLength;
        pos += commonPrefixLength;
        _range.set_begin_pos(_range.begin_pos() + commonPrefixLength);
    }

    // Added to address STREL-275
    unsigned
    getCommonPrefixLength() const
    {
        return _commonPrefixLength;
    }

    const std::vector<GermlineIndelAlleleInfo>&
    getIndelAlleles() const
    {
        return _indelAlleleInfo;
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

    /// return "totalReadDepth" summed over all samples, see method of same name in indel
    /// sample object for totalReadDepth definition.
    unsigned
    getTotalReadDepth() const
    {
        const unsigned sampleCount(getSampleCount());
        unsigned allSampleLocusDepth(0);
        for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
        {
            allSampleLocusDepth += getIndelSample(sampleIndex).getTotalReadDepth();
        }
        return allSampleLocusDepth;
    }

    /// run all internal consistency checks
    void
    assertValidity() const
    {
        const unsigned sampleCount(getSampleCount());
        assert (sampleCount == _indelSampleInfo.size());
    }

    void
    doNotGenotype()
    {
        _doNotGenotype = true;
    }

    bool
    isNotGenotyped() const
    {
        return _doNotGenotype;
    }

private:
    std::vector<GermlineIndelAlleleInfo> _indelAlleleInfo;
    std::vector<GermlineIndelSampleInfo> _indelSampleInfo;

    /// The reference range of all indel alleles at this locus:
    known_pos_range2 _range;

    unsigned _commonPrefixLength;

    /// if true genotyping is not conducted
    bool _doNotGenotype;
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
    clearEVSFeatures()
    {
        evsFeatures.clear();
        evsDevelopmentFeatures.clear();
    }

    /// production and development features used in the empirical scoring model:
    VariantScoringFeatureKeeper evsFeatures;
    VariantScoringFeatureKeeper evsDevelopmentFeatures;
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


struct GermlineSiteSampleInfo
{
    /// get the count of reads crossing the site at any mapq level
    ///
    /// TODO get a more exact definition of requrements for read to be counted here... soft-clip, overlapping frags, basecallQual = 0 , etc...
    /// TODO how does this relate to AD?
    unsigned
    getTotalReadDepth() const
    {
        return mapqTracker.count;
    }

    /// get the count of all reads crossing the site that passed mapping filters
    ///
    /// TODO get a more exact definition of requrements for read to be counted here... soft-clip, overlapping frags, basecallQual = 0 , etc...
    /// TODO how does this relate to AD?
    unsigned
    getMappedReadDepth() const
    {
        return usedBasecallCount + unusedBasecallCount;
    }

    bool
    isUsedReadCoverage() const
    {
        return (usedBasecallCount != 0);
    }

    bool
    isAnyReadCoverage() const
    {
        return (isUsedReadCoverage() or (unusedBasecallCount != 0));
    }

    void
    clear()
    {
        spanningDeletionReadCount = 0;
        usedBasecallCount = 0;
        unusedBasecallCount = 0;
        isOverlappingHomAltDeletion = false;
        mapqTracker.clear();
        ReadPosRankSum = 0;
        BaseQRankSum = 0;
        MQRankSum = 0;
        avgBaseQ = 0;
        rawPos = 0;
        strandBias = 0;

    }

    /// count of reads which have a most likely alignment containing a deletion at the site in question
    unsigned spanningDeletionReadCount = 0;

    unsigned usedBasecallCount = 0;
    unsigned unusedBasecallCount = 0;

    /// set to true when the site is overlapped by a hom deletion:
    bool isOverlappingHomAltDeletion = false;

    MapqTracker mapqTracker;

    //only meaningful for het calls
    double ReadPosRankSum = 0;  // Uses Mann-Whitney Rank Sum Test for the distance from the end of the read containing an alternate allele.
    double BaseQRankSum = 0;    // Uses Mann-Whitney Rank Sum Test for BQs (ref bases vs alternate alleles)
    double MQRankSum = 0;       // Uses Mann-Whitney Rank Sum Test for MQs (ref bases vs alternate alleles)

    double meanDistanceFromReadEdge = 0;
    double avgBaseQ = 0;
    double rawPos = 0;

    /// Sample-specific SNV bias metric:
    double strandBias = 0;
};



/// represents an site call at the level of a full VCF record, containing possibly multiple alleles
struct GermlineSiteLocusInfo : public LocusInfo
{
    typedef LocusInfo base_t;

    GermlineSiteLocusInfo(
        const unsigned sampleCount,
        const pos_t initPos,
        const uint8_t initRefBaseIndex,
        const bool initIsForcedOutput = false)
        : base_t(sampleCount, initPos),
          refBaseIndex(initRefBaseIndex),
          _siteSampleInfo(sampleCount)
    {
        isForcedOutput = initIsForcedOutput;
    }

    explicit
    GermlineSiteLocusInfo(
        const unsigned sampleCount)
        : base_t(sampleCount),
          _siteSampleInfo(sampleCount)
    {}

    bool
    isRefUnknown() const
    {
        return (refBaseIndex == BASE_ID::ANY);
    }


    void
    setSiteSampleInfo(
        const unsigned sampleIndex,
        const GermlineSiteSampleInfo& siteSampleInfo)
    {
        // ensure that no alleles are added once we start adding samples...
        _isLockAlleles = true;
        assert(sampleIndex < _siteSampleInfo.size());
        _siteSampleInfo[sampleIndex] = siteSampleInfo;
    }

    const GermlineSiteSampleInfo&
    getSiteSample(const unsigned sampleIndex) const
    {
        return _siteSampleInfo[sampleIndex];
    }

    void
    addAltSiteAllele(
        const BASE_ID::index_t baseId)
    {
        assert(not _isLockAlleles);

        _siteAlleleInfo.emplace_back(baseId);
        incrementAltAlleleCount();
    }

    std::vector<GermlineSiteAlleleInfo>&
    getSiteAlleles()
    {
        return _siteAlleleInfo;
    }

    const std::vector<GermlineSiteAlleleInfo>&
    getSiteAlleles() const
    {
        return _siteAlleleInfo;
    }

    /// return "totalReadDepth" summed over all samples, see method of same name in site
    /// sample object for totalReadDepth definition.
    unsigned
    getTotalReadDepth() const
    {
        const unsigned sampleCount(getSampleCount());
        unsigned allSampleLocusDepth(0);
        for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
        {
            allSampleLocusDepth += getSiteSample(sampleIndex).getTotalReadDepth();
        }
        return allSampleLocusDepth;
    }

    /// test whether known QUAL value should be written for this locus
    bool
    isQual() const
    {
        if (isRefUnknown()) return false;

        // test for at least one non-empty sample:
        const unsigned sampleCount(getSampleCount());
        for (unsigned sampleIndex(0); sampleIndex< sampleCount; ++sampleIndex)
        {
            const auto& sampleInfo(getSample(sampleIndex));
            const auto& siteSample(getSiteSample(sampleIndex));
            if (siteSample.isUsedReadCoverage() and (not siteSample.isOverlappingHomAltDeletion) and sampleInfo.isVariant())
            {
                return true;
            }
        }
        return false;
    }

    /// test whether known GQX value should be written for this locus/sample
    bool
    is_gqx(const unsigned sampleIndex) const
    {
        if (isRefUnknown()) return false;

        const auto& siteSample(getSiteSample(sampleIndex));
        return (siteSample.isUsedReadCoverage() and (not siteSample.isOverlappingHomAltDeletion));
    }

    void
    clear()
    {
        base_t::clear();
        refBaseIndex = BASE_ID::ANY;
        hpol = 0;
        isForcedOutput = false;
        _siteAlleleInfo.clear();

        for (auto& siteSample : _siteSampleInfo)
        {
            siteSample.clear();
        }
    }

    uint8_t refBaseIndex = BASE_ID::ANY;

    unsigned hpol = 0;

    bool isForcedOutput = false;

private:
    std::vector<GermlineSiteAlleleInfo> _siteAlleleInfo;
    std::vector<GermlineSiteSampleInfo> _siteSampleInfo;
};

std::ostream& operator<<(std::ostream& os,const GermlineSiteLocusInfo& si);


/// specify that calling model is diploid
struct GermlineDiploidSiteLocusInfo : public GermlineSiteLocusInfo
{
    GermlineDiploidSiteLocusInfo(
        const gvcf_deriv_options& gvcfDerivedOptions,
        const unsigned sampleCount,
        const pos_t init_pos,
        const uint8_t initRefBaseIndex,
        const bool is_forced_output = false)
        : GermlineSiteLocusInfo(sampleCount, init_pos, initRefBaseIndex, is_forced_output),
          evsFeatures(gvcfDerivedOptions.snvFeatureSet),
          evsDevelopmentFeatures(gvcfDerivedOptions.snvDevelopmentFeatureSet)
    {}

    GermlineDiploidSiteLocusInfo(
        const gvcf_deriv_options& gvcfDerivedOptions,
        const unsigned sampleCount)
        : GermlineSiteLocusInfo(sampleCount),
          evsFeatures(gvcfDerivedOptions.snvFeatureSet),
          evsDevelopmentFeatures(gvcfDerivedOptions.snvDevelopmentFeatureSet)
    {}

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

    void
    clearEVSFeatures()
    {
        evsFeatures.clear();
        evsDevelopmentFeatures.clear();
    }

    /// production and development features used in the empirical scoring model:
    VariantScoringFeatureKeeper evsFeatures;
    VariantScoringFeatureKeeper evsDevelopmentFeatures;
};

std::ostream& operator<<(std::ostream& os,const GermlineDiploidSiteLocusInfo& si);


/// TODO STREL-125 - transition to using regular counts structures
struct GermlineContinuousSiteSampleInfo
{
    void
    clear()
    {
        continuousTotalDepth = 0;
        continuousAlleleDepth = 0;
    }

    double
    getContinuousAlleleFrequency() const
    {
        return safeFrac(continuousAlleleDepth, continuousTotalDepth);
    }

    unsigned continuousTotalDepth = 0;
    unsigned continuousAlleleDepth = 0;
};


/// specify that variant is site/SNV and a continuous frequency calling model
struct GermlineContinuousSiteLocusInfo : public GermlineSiteLocusInfo
{
    typedef GermlineSiteLocusInfo base_t;

    GermlineContinuousSiteLocusInfo(
        const unsigned sampleCount,
        const pos_t init_pos,
        const uint8_t initRefBaseIndex,
        const bool is_forced_output = false)
        : base_t(sampleCount, init_pos, initRefBaseIndex, is_forced_output),
          _continuousSiteSampleInfo(sampleCount)
    {}

    void
    clear()
    {
        base_t::clear();

        for (auto& sample : _continuousSiteSampleInfo)
        {
            sample.clear();
        }
    }

    void
    setContinuousSiteSampleInfo(
        const unsigned sampleIndex,
        const GermlineContinuousSiteSampleInfo& siteSampleInfo)
    {
        _isLockAlleles = true;
        assert(sampleIndex < _continuousSiteSampleInfo.size());
        _continuousSiteSampleInfo[sampleIndex] = siteSampleInfo;
    }

    const GermlineContinuousSiteSampleInfo&
    getContinuousSiteSample(const unsigned sampleIndex) const
    {
        return _continuousSiteSampleInfo[sampleIndex];
    }

private:
    std::vector<GermlineContinuousSiteSampleInfo> _continuousSiteSampleInfo;
};

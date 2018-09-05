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

#include "starling_pos_processor.hh"

#include "starling_continuous_variant_caller.hh"
#include "blt_common/ref_context.hh"
#include "blt_util/log.hh"
#include "blt_util/prob_util.hh"
#include "common/Exceptions.hh"
#include "starling_common/AlleleGroupGenotype.hh"
#include "starling_common/AlleleReportInfoUtil.hh"
#include "starling_common/indel_util.hh"
#include "starling_common/OrthogonalVariantAlleleCandidateGroup.hh"
#include "starling_common/OrthogonalVariantAlleleCandidateGroupUtil.hh"
#include "gvcf_locus_info.hh"

#include <iomanip>



starling_pos_processor::
starling_pos_processor(
    const starling_options& opt,
    const starling_deriv_options& dopt,
    const reference_contig_segment& ref,
    const starling_streams& fileStreams,
    RunStatsManager& statsManager)
    : base_t(opt, dopt, ref, fileStreams, opt.alignFileOpt.alignmentFilenames.size(), statsManager),
      _opt(opt),
      _dopt(dopt),
      _streams(fileStreams)
{
    const unsigned sampleCount(getSampleCount());
    assert(_streams.getSampleNames().size() == sampleCount);

    // setup gvcf aggregator
    if (_opt.gvcf.is_gvcf_output())
    {
        _gvcfer.reset(new gvcf_aggregator(
                          _opt, _dopt, _streams, ref, _nocompress_regions, _callRegions, sampleCount));
    }

    // setup indel buffer samples:
    {
        for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
        {
            sample_info& sif(sample(sampleIndex));
            getIndelBuffer().registerSample(sif.estdepth_buff, sif.estdepth_buff_tier2, true);
        }

        getIndelBuffer().finalizeSamples();
    }
}



void
starling_pos_processor::
resetRegion(
    const std::string& chromName,
    const known_pos_range2& reportRegion)
{
    base_t::resetRegionBase(chromName, reportRegion);

    assert(_gvcfer);
    _gvcfer->resetRegion(chromName, reportRegion);

    // setup indel buffer max depth:
    {
        double maxIndelCandidateDepthSumOverNormalSamples(-1.);
        if (_dopt.gvcf.is_max_depth())
        {
            if (_opt.max_candidate_indel_depth_factor > 0.)
            {
                maxIndelCandidateDepthSumOverNormalSamples = (_opt.max_candidate_indel_depth_factor *
                                                              _gvcfer->getMaxDepth());
            }
        }

        if (_opt.max_candidate_indel_depth > 0.)
        {
            if (maxIndelCandidateDepthSumOverNormalSamples > 0.)
            {
                maxIndelCandidateDepthSumOverNormalSamples = std::min(maxIndelCandidateDepthSumOverNormalSamples,
                                                                      static_cast<double>(_opt.max_candidate_indel_depth));
            }
            else
            {
                maxIndelCandidateDepthSumOverNormalSamples = _opt.max_candidate_indel_depth;
            }
        }

        getIndelBuffer().setMaxCandidateDepth(maxIndelCandidateDepthSumOverNormalSamples);
    }
}



void
starling_pos_processor::
insert_nocompress_region(
    const known_pos_range2& range)
{
    _stagemanPtr->validate_new_pos_value(range.begin_pos(),STAGE::READ_BUFFER);
    _nocompress_regions.addRegion(range);
    _is_skip_process_pos=false;
}



void
starling_pos_processor::
reset()
{
    base_t::reset();

    if (_opt.gvcf.is_gvcf_output())
    {
        _gvcfer->reset();
    }
    _nocompress_regions.clear();
    _variantLocusAlreadyOutputToPos = -1;
}



void
starling_pos_processor::
process_pos_snp(const pos_t pos)
{
    try
    {
        const unsigned sampleCount(getSampleCount());

        const bool isForcedOutput(is_forced_output_pos(pos));

        const bool isSkippable(!isForcedOutput);

        if (isSkippable)
        {
            bool isZeroCoverage(true);
            for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
            {
                const sample_info& sif(sample(sampleIndex));
                const CleanedPileup& cpi(sif.cleanedPileup);
                const snp_pos_info& pi(cpi.rawPileup());

                if (not pi.calls.empty())
                {
                    isZeroCoverage = false;
                    break;
                }
            }

            if (isZeroCoverage) return;
        }

        // prep step 1) clean pileups in all samples:
        for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
        {
            _pileupCleaner.CleanPileupErrorProb(sample(sampleIndex).cleanedPileup);
        }

        if (_opt.is_bsnp_diploid())
        {
            process_pos_snp_digt(pos);
        }
        else
        {
            process_pos_snp_continuous(pos);
        }

    }
    catch (...)
    {
        log_os << "Exception caught in starling_pos_processor_base.process_pos_snp() while processing chromosome position: " << (pos+1) << "\n";
        throw;
    }
}



/// Setup siteSampleInfo assuming that corresponding sampleInfo has already been initialized
static
void
updateSiteSampleInfo(
    const starling_options& opt,
    const unsigned sampleIndex,
    const CleanedPileup& cpi,
    const bool isOverlappingHomAltDeletion,
    const double strandBias,
    GermlineSiteLocusInfo& locus)
{
    const snp_pos_info& good_pi(cpi.cleanedPileup());

    GermlineSiteSampleInfo siteSampleInfo;

    // note - these two values related to overlapping deletions come from diff sources, one is based on the
    // most likely genotype of overlapping indels, and the other is counts of overlapping reads.
    siteSampleInfo.isOverlappingHomAltDeletion = isOverlappingHomAltDeletion;

    siteSampleInfo.spanningDeletionReadCount = good_pi.spanningDeletionReadCount;

    siteSampleInfo.usedBasecallCount = cpi.usedBasecallCount();
    siteSampleInfo.unusedBasecallCount = cpi.unusedBasecallCount();

    // MQ is computed/reported whether EVS features are needed or not, it is also used by EVS
    const auto& pi(cpi.rawPileup());
    siteSampleInfo.mapqTracker = pi.mapqTracker;

    // Bound 'raw' strand-bias input to be less than the specified absolute value
    assert(opt.maxAbsSampleVariantStrandBias >= 0.);
    siteSampleInfo.strandBias = std::min(opt.maxAbsSampleVariantStrandBias, std::max(-opt.maxAbsSampleVariantStrandBias, strandBias));

    // add EVS feature info
    const auto& sampleInfo(locus.getSample(sampleIndex));
    if (locus.isForcedOutput or sampleInfo.isVariant())
    {
        // calculate empirical scoring metrics
        if (opt.is_compute_germline_scoring_metrics())
        {
            siteSampleInfo.ReadPosRankSum = pi.get_read_pos_ranksum();
            siteSampleInfo.MQRankSum = pi.get_mq_ranksum();
            siteSampleInfo.BaseQRankSum = pi.get_baseq_ranksum();

            siteSampleInfo.meanDistanceFromReadEdge = pi.distanceFromReadEdge.mean();
            siteSampleInfo.rawPos = pi.get_raw_pos();
            siteSampleInfo.avgBaseQ = pi.get_raw_baseQ();
        }
    }

    locus.setSiteSampleInfo(sampleIndex, siteSampleInfo);
}



static
void
computeSampleDiploidSiteGenotype(
    const starling_base_options& opt,
    const starling_deriv_options& dopt,
    const starling_pos_processor::sample_info& sif,
    const unsigned ploidy,
    diploid_genotype& dgt)
{
    const extended_pos_info& good_epi(sif.cleanedPileup.getExtendedPosInfo());
    dgt.ploidy=ploidy;
    dopt.pdcaller().position_snp_call_pprob_digt(
        opt, good_epi, dgt, opt.is_all_sites());
}



/// Translate a base index into the an allele index for the locus
///
/// site loci contain alleles in vcf print order of REF, ALT0, ALT1, etc...
/// this translates the base into the loci's allele index order
static
uint8_t
translateBaseIndexToAlleleIndex(
    const uint8_t baseIndex,
    const GermlineSiteLocusInfo& locus)
{
    const auto& siteAlleles(locus.getSiteAlleles());
    uint8_t alleleIndex(0);
    if (locus.refBaseIndex == baseIndex)
    {
        return alleleIndex;
    }
    alleleIndex++;

    for (const auto& allele : siteAlleles)
    {
        if (allele.baseIndex == baseIndex) break;
        alleleIndex++;
    }

    if (alleleIndex > siteAlleles.size())
    {
        using namespace illumina::common;
        std::ostringstream oss;
        oss << "Can't find allele " << id_to_base(baseIndex) << " expected to be present in site locus";
        BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
    }

    return alleleIndex;
}



/// Translate older starling "4-allele" genotype index to current VcfGenotype associated
/// with the specified site locus
static
void
translateDigtToVcfGenotype(
    const int ploidy,
    const unsigned digtIndex,
    const GermlineSiteLocusInfo& locus,
    VcfGenotype& vcfGt)
{
    if (ploidy == 1)
    {
        const uint8_t allele0Index(translateBaseIndexToAlleleIndex(DIGT::get_allele(digtIndex, 0), locus));
        vcfGt.setGenotypeFromAlleleIndices(allele0Index);
    }
    else if (ploidy == 2)
    {
        const uint8_t allele0Index(translateBaseIndexToAlleleIndex(DIGT::get_allele(digtIndex, 0), locus));
        const uint8_t allele1Index(translateBaseIndexToAlleleIndex(DIGT::get_allele(digtIndex, 1), locus));
        vcfGt.setGenotypeFromAlleleIndices(allele0Index, allele1Index);
    }
    else
    {
        assert(false and "Unexpected ploidy");
    }
}



static
void
updateSnvLocusWithSampleInfo(
    const starling_options& opt,
    const starling_pos_processor::sample_info& sif,
    const unsigned callerPloidy,
    const unsigned groupLocusPloidy,
    const diploid_genotype& dgt,
    const unsigned sampleIndex,
    const ActiveRegionId activeRegionId,
    const CandidateSnvBuffer& candidateSnvBuffer,
    GermlineDiploidSiteLocusInfo& locus,
    double& homRefLogProb)
{
    auto& sampleInfo(locus.getSample(sampleIndex));
    sampleInfo.setPloidy(callerPloidy);

    sampleInfo.setActiveRegionId(activeRegionId);

    bool isOverlappingHomAltDeletion(false);
    if (groupLocusPloidy == 0)
    {
        const int locusPloidyAdjustment(sif.cleanedPileup.rawPileup().spanningIndelPloidyModification);
        isOverlappingHomAltDeletion=(locusPloidyAdjustment < 0);
    }

    const CleanedPileup& cpi(sif.cleanedPileup);

    if (cpi.usedBasecallCount() != 0)
    {
        // the principle of this filter is that there's supposed to be no coverage here
        // we make an exception for sites inside of homalt deletions, maybe we shouldn't?
        if ((groupLocusPloidy == 0) and (not isOverlappingHomAltDeletion))
        {
            sampleInfo.setPloidyConflict();
        }
    }

    if     (locus.isRefUnknown() or (cpi.usedBasecallCount() == 0) or isOverlappingHomAltDeletion)
    {
        sampleInfo.genotypeQuality = 0;
        sampleInfo.maxGenotypeIndex.setGenotypeFromAlleleIndices();

        sampleInfo.genotypeQualityPolymorphic=0;
        sampleInfo.maxGenotypeIndexPolymorphic.setGenotypeFromAlleleIndices();
        sampleInfo.gqx=0;
    }
    else
    {
        sampleInfo.genotypeQuality = dgt.genome.max_gt_qphred;
        translateDigtToVcfGenotype(callerPloidy, dgt.genome.max_gt, locus, sampleInfo.maxGenotypeIndex);

        sampleInfo.genotypeQualityPolymorphic = dgt.poly.max_gt_qphred;
        translateDigtToVcfGenotype(callerPloidy, dgt.poly.max_gt, locus, sampleInfo.maxGenotypeIndexPolymorphic);

        sampleInfo.setGqx();
    }

    if (not locus.isRefUnknown())
    {
        // set PL values:
        const auto& siteAlleles(locus.getSiteAlleles());
        const uint8_t altAlleleCount(siteAlleles.size());
        const uint8_t fullAlleleCount(altAlleleCount + 1);
        const bool isAltAlleles(altAlleleCount > 0);
        if (isAltAlleles)
        {
            auto alleleIndexToBaseIndex = [&](const uint8_t alleleIndex)
            {
                if (alleleIndex == 0) return locus.refBaseIndex;
                return static_cast<uint8_t>(siteAlleles[alleleIndex-1].baseIndex);
            };

            // number of PL fields required:
            //const unsigned genotypeCount(VcfGenotypeUtil::getGenotypeCount(callerPloidy, fullAlleleCount));
            sampleInfo.genotypePhredLoghood.setPloidy(callerPloidy);

            if (callerPloidy == 1)
            {
                for (unsigned allele0Index(0); allele0Index < fullAlleleCount; ++allele0Index)
                {
                    const uint8_t base0Index(alleleIndexToBaseIndex(allele0Index));
                    sampleInfo.genotypePhredLoghood.getGenotypeLikelihood(allele0Index) =
                        dgt.phredLoghood[base0Index];
                }
            }
            else if (callerPloidy == 2)
            {
                for (unsigned allele1Index(0); allele1Index < fullAlleleCount; ++allele1Index)
                {
                    for (unsigned allele0Index(0); allele0Index <= allele1Index; ++allele0Index)
                    {
                        const uint8_t base0Index(alleleIndexToBaseIndex(allele0Index));
                        const uint8_t base1Index(alleleIndexToBaseIndex(allele1Index));
                        const unsigned digtGenotypeIndex(DIGT::get_gt_with_alleles(base0Index, base1Index));
                        sampleInfo.genotypePhredLoghood.getGenotypeLikelihood(allele0Index, allele1Index) =
                            dgt.phredLoghood[digtGenotypeIndex];
                    }
                }
            }
            else
            {
                assert(false and "Unexpected ploidy");
            }
        }

        // update AD counts:
        {
            uint8_t baseIndexToAlleleIndex[N_BASE];
            {
                std::fill(std::begin(baseIndexToAlleleIndex), std::end(baseIndexToAlleleIndex), fullAlleleCount);
                unsigned alleleIndex(0);
                baseIndexToAlleleIndex[locus.refBaseIndex] = alleleIndex;
                for (const auto& allele : siteAlleles)
                {
                    alleleIndex++;
                    baseIndexToAlleleIndex[allele.baseIndex] = alleleIndex;
                }
            }

            sampleInfo.supportCounts.setAltCount(altAlleleCount);

            const snp_pos_info& good_pi(cpi.cleanedPileup());
            for (const auto& call : good_pi.calls)
            {
                if (call.base_id == BASE_ID::ANY) continue;
                const uint8_t alleleIndex(baseIndexToAlleleIndex[call.base_id]);
                if (alleleIndex == fullAlleleCount) continue;
                sampleInfo.supportCounts.getCounts(call.is_fwd_strand).incrementAlleleCount(alleleIndex);
            }
        }

        // update homref prob for QUAL
        homRefLogProb += std::log(dgt.genome.ref_pprob);

        /// TODO STREL-125 find a way to restore strand bias feature
        // allele.strandBias=dgt.strand_bias;
    }

    updateSiteSampleInfo(opt, sampleIndex, cpi, isOverlappingHomAltDeletion, dgt.strand_bias, locus);

    // set complex allele ids
    auto& maxGt = sampleInfo.max_gt();
    if (maxGt.getPloidy() >= 1)
    {
        auto allele0Index(maxGt.getAllele0Index());
        if (allele0Index > 0)
        {
            maxGt.setAllele0HaplotypeId(
                candidateSnvBuffer.getHaplotypeId(sampleIndex,
                                                  locus.pos, locus.getSiteAlleles()[allele0Index-1].baseIndex)
            );
        }
    }
    if (maxGt.getPloidy() == 2)
    {
        auto allele1Index(maxGt.getAllele1Index());
        if (allele1Index > 0)
        {
            maxGt.setAllele1HaplotypeId(
                candidateSnvBuffer.getHaplotypeId(sampleIndex,
                                                  locus.pos, locus.getSiteAlleles()[allele1Index-1].baseIndex)
            );
        }
    }
    maxGt.addAltAlleleHaplotypeCountRatio(candidateSnvBuffer.getAltHaplotypeCountRatio(sampleIndex, locus.pos));
}



void
starling_pos_processor::
getSiteAltAlleles(
    const uint8_t refBaseIndex,
    const std::vector<diploid_genotype>& allDgt,
    std::vector<uint8_t>& altAlleles) const
{
    const unsigned sampleCount(getSampleCount());

    // rank alleles in each sample, sum score for each sample: first=2, second=1, not present in top PLOIDY alleles = 0
    std::array<unsigned, N_BASE> alleleRank;
    std::fill(std::begin(alleleRank), std::end(alleleRank), 0);

    std::array<double, N_BASE> sampleBaseCounts;
    for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
    {
        const auto& cpi(sample(sampleIndex).cleanedPileup);
        const auto& good_pi(cpi.cleanedPileup());
        good_pi.getBasecallCounts(sampleBaseCounts);

        static const double minAlleleFraction(0.10);
        unsigned minCount(0);
        {
            for (unsigned baseIndex(0); baseIndex < N_BASE; ++baseIndex)
            {
                minCount += sampleBaseCounts[baseIndex];
            }
            minCount *= minAlleleFraction;
            minCount = std::max(1u, minCount);
        }

        const auto ploidy(allDgt[sampleIndex].ploidy);
        for (int ploidyIndex(0); ploidyIndex < ploidy; ++ploidyIndex)
        {
            unsigned maxBaseIndex(0);
            for (unsigned baseIndex(1); baseIndex < N_BASE; ++baseIndex)
            {
                if (sampleBaseCounts[baseIndex] > sampleBaseCounts[maxBaseIndex])
                {
                    maxBaseIndex = baseIndex;
                }
            }
            if (sampleBaseCounts[maxBaseIndex] >= minCount)
            {
                alleleRank[maxBaseIndex] += (2 - ploidyIndex);
            }
            sampleBaseCounts[maxBaseIndex] = 0;
        }
    }

    // translate rank score into order altAllele list:
    bool isBaseAdded[N_BASE];
    for (unsigned baseIndex(0); baseIndex < N_BASE; ++baseIndex)
    {
        isBaseAdded[baseIndex] = (baseIndex == refBaseIndex);
    }
    while (true)
    {
        unsigned maxBaseIndex(0);
        for (unsigned baseIndex(1); baseIndex < N_BASE; ++baseIndex)
        {
            if (alleleRank[baseIndex] > alleleRank[maxBaseIndex])
            {
                maxBaseIndex = baseIndex;
            }
        }
        if (alleleRank[maxBaseIndex] == 0) break;
        alleleRank[maxBaseIndex] = 0;
        if (maxBaseIndex == refBaseIndex) continue;
        altAlleles.push_back(maxBaseIndex);
        isBaseAdded[maxBaseIndex] = true;
    }

    // go through most likely genotypes and add any remaining bases in essentially random order
    for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
    {
        const auto& dgt(allDgt[sampleIndex]);
        const int ploidy(dgt.ploidy);

        auto checkGtChrom = [&](const unsigned genotypeIndex, const unsigned chromIndex)
        {
            const uint8_t baseIndex(DIGT::get_allele(genotypeIndex, chromIndex));
            if (isBaseAdded[baseIndex]) return;
            altAlleles.push_back(baseIndex);
            isBaseAdded[baseIndex] = true;
        };

        auto checkGt = [&](const unsigned genotypeIndex)
        {
            if (ploidy == 1)
            {
                checkGtChrom(genotypeIndex, 0);
            }
            else if (ploidy == 2)
            {
                checkGtChrom(genotypeIndex, 0);
                checkGtChrom(genotypeIndex, 1);
            }
            else
            {
                assert(false and "Unexpected ploidy");
            }
        };

        checkGt(dgt.poly.max_gt);
        checkGt(dgt.genome.max_gt);
    }
}



void
starling_pos_processor::
process_pos_snp_digt(
    const pos_t pos)
{
    const unsigned sampleCount(getSampleCount());
    const bool isForcedOutput(is_forced_output_pos(pos));

    /// TODO currently using old "4-allele" genotyping system and then reducing down to number
    /// of called alleles, transition this to work more like current indel model, where up to
    /// 4 alleles are nominated as candidates and genotyping is based on the candidate alleles
    /// only.

    // -----------------------------------------------
    // prep data for site locus creation:
    //

    // prep step 2) setup ploidy arrays:
    std::vector<int> groupLocusPloidy;
    std::vector<int> callerPloidy;
    for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
    {
        // groupLocusPloidy of 0 is treated as a special case, if this happens the
        // entire calling method reverts to a ploidy of 2 for the sample, but the
        // locus ploidy is passed into the gVCF writer as 0. The gVCF writer can
        // decide what to do with this information from there.
        //
        const int regionPloidy(get_ploidy(pos, sampleIndex));
        const int locusPloidyAdjustment(sample(sampleIndex).cleanedPileup.rawPileup().spanningIndelPloidyModification);
        const int ploidy = std::max(0, regionPloidy+locusPloidyAdjustment);
        groupLocusPloidy.push_back(ploidy);
        callerPloidy.push_back((ploidy == 0) ? 2 : ploidy);
    }

    // prep step 3) compute diploid genotype object using older "4-allele" model:
    std::vector<diploid_genotype> allDgt(sampleCount);
    for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
    {
        computeSampleDiploidSiteGenotype(
            _opt, _dopt, sample(sampleIndex), callerPloidy[sampleIndex], allDgt[sampleIndex]);
    }

    // prep step 4) rank each allele in each sample, allowing up to ploidy alleles.
    //              approximate an aggregate rank over all samples:
    const uint8_t refBaseIndex(base_to_id(_ref.get_base(pos)));
    std::vector<uint8_t> altAlleles;
    if (refBaseIndex != BASE_ID::ANY)
    {
        getSiteAltAlleles(refBaseIndex, allDgt, altAlleles);
    }

    // -----------------------------------------------
    // create site locus object:
    //
    std::unique_ptr<GermlineDiploidSiteLocusInfo> locusPtr(new GermlineDiploidSiteLocusInfo(_dopt.gvcf, sampleCount, pos, refBaseIndex, isForcedOutput));


    // add all candidate alternate alleles:
    for (const auto baseId : altAlleles)
    {
        locusPtr->addAltSiteAllele(static_cast<BASE_ID::index_t>(baseId));
    }

    double homRefLogProb(0);
    ActiveRegionId activeRegionId(getActiveRegionDetector().getActiveRegionId(pos));
    for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
    {
        updateSnvLocusWithSampleInfo(
            _opt, sample(sampleIndex), callerPloidy[sampleIndex], groupLocusPloidy[sampleIndex],
            allDgt[sampleIndex], sampleIndex, activeRegionId, getCandidateSnvBuffer(), *locusPtr, homRefLogProb);
    }

    // add sample-independent info:
    locusPtr->anyVariantAlleleQuality = ln_error_prob_to_qphred(homRefLogProb);

    if (locusPtr->isVariantLocus())
    {
        // hpol filter
        locusPtr->hpol = get_snp_hpol_size(pos,_ref);
    }

    //Add site to gvcf
    _gvcfer->add_site(std::move(locusPtr));
}



static
void
updateContinuousSiteSampleInfo(
    const unsigned sampleIndex,
    const unsigned alleleIndex, ///< TODO STREL-125 TMP!!!!!!
    GermlineContinuousSiteLocusInfo& locus)
{
    const auto& sampleInfo(locus.getSample(sampleIndex));
    const auto& siteSampleInfo(locus.getSiteSample(sampleIndex));

    GermlineContinuousSiteSampleInfo siteContinuousSampleInfo;

    siteContinuousSampleInfo.continuousTotalDepth = siteSampleInfo.spanningDeletionReadCount;
    const auto& fwd(sampleInfo.supportCounts.getCounts(true));
    const auto& rev(sampleInfo.supportCounts.getCounts(false));
    siteContinuousSampleInfo.continuousTotalDepth += fwd.totalConfidentCounts() + rev.totalConfidentCounts();
    siteContinuousSampleInfo.continuousTotalDepth += fwd.nonConfidentCount + rev.nonConfidentCount;
    siteContinuousSampleInfo.continuousAlleleDepth = fwd.confidentAlleleCount(alleleIndex) + rev.confidentAlleleCount(alleleIndex);

    locus.setContinuousSiteSampleInfo(sampleIndex, siteContinuousSampleInfo);
}



/// strand bias values summed over all samples
struct StrandBiasCounts
{
    void
    merge(const StrandBiasCounts& rhs)
    {
        fwdAlt += rhs.fwdAlt;
        revAlt += rhs.revAlt;
        fwdOther += rhs.fwdOther;
        revOther += rhs.revOther;
    }

    unsigned fwdAlt = 0;
    unsigned revAlt = 0;
    unsigned fwdOther = 0;
    unsigned revOther = 0;
};



static
void
updateContinuousSnvLocusWithSampleInfo(
    const starling_options& opt,
    starling_pos_processor::sample_info& sif,
    const unsigned sampleIndex,
    std::vector<StrandBiasCounts>& strandBiasCounts,
    GermlineContinuousSiteLocusInfo& locus)
{
    const auto& siteAlleles(locus.getSiteAlleles());
    const bool isRefAllele(siteAlleles.empty());

    const uint8_t altAlleleCount(siteAlleles.size());
    const uint8_t fullAlleleCount(altAlleleCount+1);

    const CleanedPileup& cpi(sif.cleanedPileup);
    const snp_pos_info& good_pi(cpi.cleanedPileup());

    auto& sampleInfo(locus.getSample(sampleIndex));

    // update AD counts:
    {
        uint8_t baseIndexToAlleleIndex[N_BASE];
        {
            std::fill(std::begin(baseIndexToAlleleIndex), std::end(baseIndexToAlleleIndex), fullAlleleCount);
            unsigned alleleIndex(0);
            baseIndexToAlleleIndex[locus.refBaseIndex] = alleleIndex;
            for (const auto& allele : siteAlleles)
            {
                alleleIndex++;
                baseIndexToAlleleIndex[allele.baseIndex] = alleleIndex;
            }
        }

        sampleInfo.supportCounts.setAltCount(altAlleleCount);
        for (const auto& call : good_pi.calls)
        {
            if (call.base_id==BASE_ID::ANY) continue;
            const uint8_t alleleIndex(baseIndexToAlleleIndex[call.base_id]);
            auto& strandCounts(sampleInfo.supportCounts.getCounts(call.is_fwd_strand));
            if (alleleIndex==fullAlleleCount)
            {
                strandCounts.nonConfidentCount++;
            }
            else
            {
                strandCounts.incrementAlleleCount(alleleIndex);
            }
        }
    }

    // update sample strand bias and locus-level strand bias intermediates:
    double strandBias(0.);
    if (not isRefAllele)
    {
        // get "primaryAltAllele" for the purpose of computing strand bias, for now this is fixed
        // to the first alt:
        unsigned primaryAltAlleleIndex(0);
        for (unsigned altAlleleIndex(0); altAlleleIndex < altAlleleCount; ++altAlleleIndex)
        {
            StrandBiasCounts sampleStrandBias;
            const auto& allele(siteAlleles[altAlleleIndex]);

            for (const base_call& bc : good_pi.calls)
            {
                if (bc.is_fwd_strand)
                {
                    if (bc.base_id == allele.baseIndex)
                        sampleStrandBias.fwdAlt++;
                    else
                        sampleStrandBias.fwdOther++;
                }
                else if (bc.base_id == allele.baseIndex)
                    sampleStrandBias.revAlt++;
                else
                    sampleStrandBias.revOther++;
            }

            if (altAlleleIndex == primaryAltAlleleIndex)
            {
                strandBias = starling_continuous_variant_caller::strandBias(sampleStrandBias.fwdAlt,
                                                                            sampleStrandBias.revAlt,
                                                                            sampleStrandBias.fwdOther,
                                                                            sampleStrandBias.revOther);
            }

            auto& sbcounts(strandBiasCounts[altAlleleIndex]);
            sbcounts.merge(sampleStrandBias);
        }
    }

    static const bool isOverlappingHomAltDeletion(false);
    updateSiteSampleInfo(opt, sampleIndex, cpi, isOverlappingHomAltDeletion, strandBias, locus);

    static const uint8_t variantAlleleIndex(1);
    const uint8_t updateAlleleIndex(isRefAllele ? 0 : variantAlleleIndex);
    updateContinuousSiteSampleInfo(sampleIndex, updateAlleleIndex, locus);

    {
        const auto& continuousSiteSampleInfo(locus.getContinuousSiteSample(sampleIndex));
        const double alleleFrequency(continuousSiteSampleInfo.getContinuousAlleleFrequency());

        {
            // use diploid gt codes as a convenient way to summarize the continuous variant calls:
            static const VcfGenotype homrefGtIndex(0, 0);
            static const VcfGenotype hetGtIndex(0, 1);
            static const VcfGenotype homGtIndex(1, 1);

            auto getGtIndex = [&]() -> VcfGenotype
            {
                if (isRefAllele)
                    return homrefGtIndex;
                else if (alleleFrequency >= (1. - opt.min_het_vf))
                    return homGtIndex;
                else if (alleleFrequency < opt.min_het_vf)
                    return homrefGtIndex; // STAR-66 - desired behavior
                else
                    return hetGtIndex;
            };

            sampleInfo.maxGenotypeIndexPolymorphic = getGtIndex();
        }

        sampleInfo.gqx = sampleInfo.genotypeQualityPolymorphic =
                             starling_continuous_variant_caller::getAlleleSequencingErrorQscore(
                                 continuousSiteSampleInfo.continuousAlleleDepth,
                                 continuousSiteSampleInfo.continuousTotalDepth,
                                 opt.continuousFrequencyCallerExpectedObservationQuality,
                                 opt.continuousFrequencyCallerMaxQscore);
    }
}



static
void
updateContinuousSnvLocusInfo(
    const starling_options& /*opt*/,
    const reference_contig_segment& ref,
    const pos_t pos,
    const std::vector<StrandBiasCounts>& strandBiasCounts,
    GermlineContinuousSiteLocusInfo& locus)
{
    const unsigned sampleCount(locus.getSampleCount());

    auto& siteAlleles(locus.getSiteAlleles());
    const unsigned alleleCount(siteAlleles.size());

    // set sample-independent info:
    locus.hpol = get_snp_hpol_size(pos, ref);

    // update locus-level strand bias for each allele:
    for (unsigned alleleIndex(0); alleleIndex < alleleCount; ++alleleIndex)
    {
        auto& allele(siteAlleles[alleleIndex]);
        const auto& sbcounts(strandBiasCounts[alleleIndex]);
        allele.strandBias = starling_continuous_variant_caller::strandBias(sbcounts.fwdAlt, sbcounts.revAlt,
                                                                           sbcounts.fwdOther, sbcounts.revOther);
    }

    // get the qual score:
    locus.anyVariantAlleleQuality = 0;
    if (alleleCount != 0)
    {
        for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
        {
            const auto& sampleInfo(locus.getSample(sampleIndex));
            locus.anyVariantAlleleQuality = std::max(locus.anyVariantAlleleQuality,
                                                     sampleInfo.genotypeQualityPolymorphic);
        }
    }
}



void
starling_pos_processor::
process_pos_snp_continuous(const pos_t pos)
{
    const unsigned sampleCount(getSampleCount());
    const bool isForcedOutput(is_forced_output_pos(pos));

    const uint8_t refBaseIndex(base_to_id(_ref.get_base(pos)));

    // report one locus (ie. vcf record) per alt allele in continuous mode
    bool isAnySiteOutputAtPosition(false);

    auto addBase = [&](const uint8_t baseIndex, const bool isForcedOutputUsed)
    {
        const bool isRefAllele(baseIndex == refBaseIndex);
        std::unique_ptr<GermlineContinuousSiteLocusInfo> locusPtr(new GermlineContinuousSiteLocusInfo(
                                                                      sampleCount, pos, refBaseIndex, isForcedOutput));

        // setup alt allele first:
        if (not isRefAllele)
        {
            locusPtr->addAltSiteAllele(static_cast<BASE_ID::index_t>(baseIndex));
        }

        // set some sample-dependent info
        const auto& siteAlleles(locusPtr->getSiteAlleles());
        const auto alleleCount(siteAlleles.size());
        std::vector<StrandBiasCounts> strandBiasCounts(alleleCount);

        for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
        {
            updateContinuousSnvLocusWithSampleInfo(
                _opt, sample(sampleIndex), sampleIndex, strandBiasCounts, *locusPtr);
        }

        // determine if this locus (in continuous case locus TEMPORARILY means "allele") is printable:
        bool isReportableLocus(isForcedOutputUsed);
        if ((not isReportableLocus) and (not isRefAllele))
        {
            for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
            {
                const auto& continuousSiteSampleInfo(locusPtr->getContinuousSiteSample(sampleIndex));
                const double alleleFrequency(continuousSiteSampleInfo.getContinuousAlleleFrequency());

                if (alleleFrequency > _opt.min_het_vf)
                {
                    isReportableLocus = true;
                    break;
                }
            }
        }
        if (not isReportableLocus) return;

        // set sample-independent info:
        updateContinuousSnvLocusInfo(
            _opt, _ref, pos, strandBiasCounts,*locusPtr);

        isAnySiteOutputAtPosition = true;
        _gvcfer->add_site(std::move(locusPtr));
    };

    for (unsigned baseIndex(0); baseIndex < N_BASE; ++baseIndex)
    {
        if (baseIndex == refBaseIndex) continue;
        addBase(baseIndex,isForcedOutput);
    }

    // ensure that at least one base is added for site
    if (not isAnySiteOutputAtPosition)
    {
        addBase(refBaseIndex,true);
    }
}



void
starling_pos_processor::
process_pos_indel(
    const pos_t pos,
    const bool isPosPrecedingReportableRange)
{
    if (_opt.is_bsnp_diploid())
    {
        // STREL-392 for this variant type to synchronize correctly across independent processes setup to call adjacent
        // regions of the genome, we need to go through a dummy calling operation for approx 'maxIndelSize' bases
        // before starting to report calls, so this methods is called even when "isPosPrecedingReportableRange"
        // is true. This helps us to correctly handle groups of large/overlapping deletions which might span the
        // segmented genome boundary between two processes. Without this step is is possible to duplicate or drop
        // indel calls across a segment boundary.
        process_pos_indel_digt(pos);
    }
    else
    {
        if (not isPosPrecedingReportableRange)
        {
            process_pos_indel_continuous(pos);
        }
    }
}



static
void
addIndelAllelesToLocus(
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    const bool isForcedOutput,
    GermlineIndelLocusInfo& locus)
{
    const unsigned nonRefAlleleCount(alleleGroup.size());
    for (unsigned nonRefAlleleIndex(0); nonRefAlleleIndex < nonRefAlleleCount; ++nonRefAlleleIndex)
    {
        const IndelKey& indelKey(alleleGroup.key(nonRefAlleleIndex));
        const IndelData& indelData(alleleGroup.data(nonRefAlleleIndex));

        // right now, a locus-level forced output flag should only correspond to forced alleles:
        assert ((not isForcedOutput) || indelData.isForcedOutput);

        locus.addAltIndelAllele(indelKey, indelData);
    }
}



/// Get the length of the common prefix of the reference and all ALT alleles
///
/// STREL-275: Added to address the bug that the REF and ALTs share a common prefix of length >1.
/// This problem happened when an internal and external deletions overlap and
/// the flanking sequences happened to be the same. E.g.
///
/// AATATATT (REF)
/// A----ATT (Internal deletion)
/// AATA---- (External deletion)
///
/// In this case, the vcf record becomes unnormalized, like
///
/// chr20 pos . AATATATT AATT,AATA ...
///
/// whereas the correct record should be
///
/// chr20 pos+2 . TATATT TT,TA ...
///
/// STREL-275 implements a suboptimal fix for NS5 by
/// (1) detecting and marking the locus with common prefixes and
/// (2) changing the output to remove the common prefix.
///
/// TODO: A more natural solution can be implemented later when the haplotyping is integrated up to the genotyping stage.
///
/// \param locus locus to be investigated
/// \param ref reference
/// \return the length of the common prefix of the reference and all ALT alleles
static
unsigned
getCommonPrefixLength(
    const GermlineIndelLocusInfo& locus,
    const reference_contig_segment& ref)
{
    const auto& locusRange(locus.range());

    unsigned minLenCommonPrefix(locusRange.size());
    for (const auto& indelAllele : locus.getIndelAlleles())
    {
        const auto& indelKey = indelAllele.indelKey;

        unsigned lenCommonPrefix(indelKey.pos - locusRange.begin_pos());
        bool isOutsideOfCommonPrefix = false;
        for (unsigned i(0); i<indelKey.insert_length(); ++i)
        {
            if (ref.get_base(locusRange.begin_pos()+lenCommonPrefix) != indelKey.insert_seq()[i])
            {
                isOutsideOfCommonPrefix = true;
                break;
            }
            ++lenCommonPrefix;
        }

        if (not isOutsideOfCommonPrefix)
        {
            for (pos_t pos(indelKey.right_pos()); pos<locusRange.end_pos(); ++pos)
            {
                if (ref.get_base(locusRange.begin_pos()+lenCommonPrefix) != ref.get_base(pos))
                {
                    break;
                }
                ++lenCommonPrefix;
            }
        }
        if (lenCommonPrefix < minLenCommonPrefix)
        {
            minLenCommonPrefix = lenCommonPrefix;
        }
    }

    return minLenCommonPrefix;
}



/// Get various indel stats from the pileup
static
void
addIndelSamplePileupInfo(
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    const pos_basecall_buffer& basecallBuffer,
    const GermlineIndelLocusInfo& locus,
    GermlineIndelSampleInfo& indelSampleInfo)
{
    const auto& range(locus.range());
    pos_t pileupPos(range.begin_pos()-1);
    const IndelKey& indelKey0(alleleGroup.key(0));
    if (indelKey0.type == INDEL::BP_RIGHT) pileupPos=range.end_pos();
    const snp_pos_info& spi(basecallBuffer.get_pos(pileupPos));
    indelSampleInfo.tier1Depth=spi.calls.size();
    indelSampleInfo.mapqTracker=spi.mapqTracker;
}



/// Add misc sample info from legacy sample indel report
static
void
addIndelSampleLegacyInfo(
    const starling_base_options& opt,
    const starling_deriv_options& dopt,
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    const unsigned sampleIndex,
    const pos_basecall_buffer& basecallBuffer,
    const bool isUseAltIndel,
    GermlineIndelSampleInfo& indelSampleInfo)
{
    static const bool isUseTier2Data(false);

    /// TODO STREL-125 legacy structure assumes single indel allele, get rid of this....
    const IndelKey& indelKey(alleleGroup.key(0));
    const IndelData& indelData(alleleGroup.data(0));
    const IndelSampleData& indelSampleData(indelData.getSampleData(sampleIndex));

    getAlleleSampleReportInfo(opt, dopt, indelKey, indelSampleData, basecallBuffer, isUseTier2Data,
                              isUseAltIndel, indelSampleInfo.legacyReportInfo);
}



/// Populate the indelSampleInfo object and add it to locus
///
/// This covers the update for both diploid/haploid and continuous calling models
///
/// \param[in] isUseAltIndel If true, account for indels other than indelKey while computing indel posteriors
static
void
updateIndelSampleInfo(
    const starling_base_options& opt,
    const starling_deriv_options& dopt,
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    const unsigned sampleIndex,
    const pos_basecall_buffer& basecallBuffer,
    const bool isUseAltIndel,
    GermlineIndelLocusInfo& locus)
{
    GermlineIndelSampleInfo indelSampleInfo;

    addIndelSamplePileupInfo(alleleGroup, basecallBuffer, locus, indelSampleInfo);
    addIndelSampleLegacyInfo(opt, dopt, alleleGroup, sampleIndex, basecallBuffer, isUseAltIndel, indelSampleInfo);

    locus.setIndelSampleInfo(sampleIndex, indelSampleInfo);
}



/// Make all sample-specific changes to the indel locus for the diploid/haploid calling case.
///
/// This will setup indelSampleInfo assuming that corresponding sampleInfo has already been initialized,
/// and update ploidy values in basecall buffer
///
/// \param[in,out] basecallBuffer Read previous base depth and mqpq values from basecallbuffer, and also
///                               set ploidy modifications on any bases which the indel overlaps
static
void
updateDiploidIndelLocusWithSampleInfo(
    const starling_base_options& opt,
    const starling_deriv_options& dopt,
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    const unsigned sampleIndex,
    pos_basecall_buffer& basecallBuffer,
    GermlineIndelLocusInfo& locus)
{
    const LocusSampleInfo& sampleInfo(locus.getSample(sampleIndex));
    const unsigned callerPloidy(sampleInfo.getPloidy().getPloidy());

    // set site ploidy modifications:
    auto updateSitePloidyForAlleleIndex = [&](const uint8_t alleleIndex)
    {
        if (alleleIndex == 0) return;

        const IndelKey& indelKey2(locus.getIndelAlleles()[alleleIndex - 1].indelKey);

        const pos_t beginPos(indelKey2.pos), endPos(indelKey2.right_pos());
        for (pos_t pos(beginPos); pos < endPos; ++pos)
        {
            basecallBuffer.decrementSpanningIndelPloidy(pos);
        }
    };

    const auto& maxGt(sampleInfo.max_gt());
    if (callerPloidy == 2)
    {
        updateSitePloidyForAlleleIndex(maxGt.getAllele0Index());
        updateSitePloidyForAlleleIndex(maxGt.getAllele1Index());
    }
    else if (callerPloidy == 1)
    {
        updateSitePloidyForAlleleIndex(maxGt.getAllele0Index());
    }
    else
    {
        assert(false and "Unexpected ploidy value");
    }

    static const bool is_use_alt_indel(false);
    updateIndelSampleInfo(opt,dopt, alleleGroup, sampleIndex, basecallBuffer, is_use_alt_indel, locus);
}



/// ploidy infered from argument count:
static
unsigned
getPriorIndex(
    const unsigned topAlleleIndexInSample,
    const unsigned allele0Index)
{
    using namespace AG_GENOTYPE;
    if (allele0Index == 0)
    {
        return HOMREF;
    }
    else
    {
        if (allele0Index == (topAlleleIndexInSample+1))
        {
            return HOM0;
        }
        else
        {
            return HOM1;
        }
    }
}



/// ploidy infered from argument count:
static
unsigned
getPriorIndex(
    const unsigned topAlleleIndexInSample,
    const unsigned allele0Index,
    const unsigned allele1Index)
{
    using namespace AG_GENOTYPE;
    if ((allele0Index == 0) and (allele1Index == 0))
    {
        return HOMREF;
    }
    else
    {
        if (allele0Index == allele1Index)
        {
            if (allele0Index == (topAlleleIndexInSample+1))
            {
                return HOM0;
            }
            else
            {
                return HOM1;
            }
        }
        else
        {
            if (allele0Index==0)
            {
                if (allele1Index == (topAlleleIndexInSample+1))
                {
                    return HET0;
                }
                else
                {
                    return HET1;
                }
            }
            else
            {
                return HET01;
            }
        }
    }
}



/// fill in all sample-specific locus info
///
/// includes:
/// -- compute GT lhoods, use this to produce GQ, PL, GQX, etc..
/// -- also use lhoods to produce support counts (AD/ADR/ADR)
/// -- other misc locus per-sample reporting requirements
///
/// \param contrastGroup these are alleles meant to be grouped together into an "other" category, such as that reported
///                      as the <*> allele in VCF
/// \param[in,out] basecallBuffer see updateIndelSampleInfo
static
void
updateIndelLocusWithSampleInfo(
    const starling_base_options& opt,
    const starling_deriv_options& dopt,
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    const unsigned topAlleleIndexInSample,
    const OrthogonalVariantAlleleCandidateGroup& contrastGroup,
    const starling_sample_options& sample_opt,
    const unsigned callerPloidy,
    const unsigned groupLocusPloidy,
    const unsigned sampleIndex,
    const ActiveRegionId activeRegionId,
    pos_basecall_buffer& basecallBuffer,
    GermlineIndelLocusInfo& locus,
    double& homRefLogProb)
{
    auto& sampleInfo(locus.getSample(sampleIndex));
    sampleInfo.setPloidy(callerPloidy);
    if (callerPloidy != groupLocusPloidy)
    {
        sampleInfo.setPloidyConflict();
    }

    sampleInfo.setActiveRegionId(activeRegionId);

    // Check if all nonref alleles are flagged as doNotGenotype
    // in such cases, this locus is marked doNotGenotype
    //
    // Note that do-not-genotype alleles should always appear on their own locus.
    //
    const unsigned nonrefAlleleGroupSize(alleleGroup.size());
    if ((nonrefAlleleGroupSize == 1) && (alleleGroup.data(0).doNotGenotype))
    {
        locus.doNotGenotype();
        return;
    }

    for (unsigned nonrefAlleleIndex(0); nonrefAlleleIndex<nonrefAlleleGroupSize; ++nonrefAlleleIndex)
    {
        const auto& indelData(alleleGroup.data(nonrefAlleleIndex));
        if (indelData.doNotGenotype)
        {
            assert(false && "Unexpected allele with doNotGenotype status");
        }
    }

    // score all possible genotypes from topVariantAlleleGroup for this sample
    std::vector<double> genotypeLogLhood;
    getVariantAlleleGroupGenotypeLhoodsForSample(
        opt, dopt, sample_opt, callerPloidy, sampleIndex, alleleGroup, contrastGroup,
        genotypeLogLhood, sampleInfo.supportCounts);

    // set phredLoghood:
    {
        const unsigned gtCount(genotypeLogLhood.size());
        unsigned maxGtIndex(0);
        for (unsigned gtIndex(1); gtIndex<gtCount; ++gtIndex)
        {
            if (genotypeLogLhood[gtIndex] > genotypeLogLhood[maxGtIndex]) maxGtIndex = gtIndex;
        }

        auto& samplePL(sampleInfo.genotypePhredLoghood.getGenotypeLikelihood());
        samplePL.resize(gtCount);
        for (unsigned gtIndex(0); gtIndex<gtCount; ++gtIndex)
        {
            samplePL[gtIndex] = ln_error_prob_to_qphred(genotypeLogLhood[gtIndex]-genotypeLogLhood[maxGtIndex]);
        }
    }

    //------------------------------------------------
    // compute posteriors/qualities:

    // get patternRepeatCount
    // if more than one alt allele then base this off of the most likely allele per sample
    const IndelData& allele0Data(alleleGroup.data(topAlleleIndexInSample));
    const AlleleReportInfo& indelReportInfo(allele0Data.getReportInfo());
    const unsigned repeatingPatternSize(std::max(1u,indelReportInfo.repeatUnitLength));
    const unsigned patternRepeatCount(std::max(1u,indelReportInfo.refRepeatCount));

    const ContextGenotypePriors& genotypePriors(dopt.getIndelGenotypePriors().getContextSpecificPriorSet(repeatingPatternSize, patternRepeatCount));

    const uint8_t nonRefAlleleCount(alleleGroup.size());
    const uint8_t fullAlleleCount(nonRefAlleleCount+1);

    // get polymorphic posterior
    std::vector<double> genotypePosterior(genotypeLogLhood.size());
    if (callerPloidy == 1)
    {
        static const bool isHaploid(true);
        const double* genotypeLogPrior(genotypePriors.getNAllelePolymorphic(isHaploid));
        for (unsigned allele0Index(0); allele0Index < fullAlleleCount; ++allele0Index)
        {
            const unsigned genotypeIndex(VcfGenotypeUtil::getGenotypeIndex(allele0Index));
            const unsigned priorIndex(getPriorIndex(topAlleleIndexInSample, allele0Index));
            genotypePosterior[genotypeIndex] = genotypeLogLhood[genotypeIndex] + genotypeLogPrior[priorIndex];
        }
    }
    else if (callerPloidy == 2)
    {
        static const bool isHaploid(false);
        const double* genotypeLogPrior(genotypePriors.getNAllelePolymorphic(isHaploid));
        for (unsigned allele1Index(0); allele1Index < fullAlleleCount; ++allele1Index)
        {
            for (unsigned allele0Index(0); allele0Index <= allele1Index; ++allele0Index)
            {
                const unsigned genotypeIndex(VcfGenotypeUtil::getGenotypeIndex(allele0Index, allele1Index));
                const unsigned priorIndex(getPriorIndex(topAlleleIndexInSample, allele0Index, allele1Index));
                genotypePosterior[genotypeIndex] = genotypeLogLhood[genotypeIndex] + genotypeLogPrior[priorIndex];
            }
        }
    }
    else
    {
        assert(false and "Unexpected ploidy value");
    }

    unsigned maxGenotypeIndex(0);
    normalizeLogDistro(std::begin(genotypePosterior), std::end(genotypePosterior), maxGenotypeIndex);

    sampleInfo.genotypeQualityPolymorphic = error_prob_to_qphred(prob_comp(std::begin(genotypePosterior), std::end(genotypePosterior), maxGenotypeIndex));

    setVcfGenotypeFromGenotypeIndex(callerPloidy, maxGenotypeIndex, sampleInfo.maxGenotypeIndexPolymorphic);

    // get genome posterior
    if (callerPloidy == 1)
    {
        static const bool isHaploid(true);
        const double* genotypeLogPrior(genotypePriors.getNAllele(isHaploid));
        for (unsigned allele0Index(0); allele0Index < fullAlleleCount; ++allele0Index)
        {
            const unsigned genotypeIndex(VcfGenotypeUtil::getGenotypeIndex(allele0Index));
            const unsigned priorIndex(getPriorIndex(topAlleleIndexInSample, allele0Index));
            genotypePosterior[genotypeIndex] = genotypeLogLhood[genotypeIndex] + genotypeLogPrior[priorIndex];
        }
    }
    else if (callerPloidy == 2)
    {
        static const bool isHaploid(false);
        const double* genotypeLogPrior(genotypePriors.getNAllele(isHaploid));
        for (unsigned allele1Index(0); allele1Index < fullAlleleCount; ++allele1Index)
        {
            for (unsigned allele0Index(0); allele0Index <= allele1Index; ++allele0Index)
            {
                const unsigned genotypeIndex(VcfGenotypeUtil::getGenotypeIndex(allele0Index, allele1Index));
                const unsigned priorIndex(getPriorIndex(topAlleleIndexInSample, allele0Index, allele1Index));
                genotypePosterior[genotypeIndex] = genotypeLogLhood[genotypeIndex] + genotypeLogPrior[priorIndex];
            }
        }
    }
    else
    {
        assert(false and "Unexpected ploidy value");
    }

    normalizeLogDistro(std::begin(genotypePosterior), std::end(genotypePosterior), maxGenotypeIndex);

    sampleInfo.genotypeQuality = error_prob_to_qphred(prob_comp(std::begin(genotypePosterior), std::end(genotypePosterior), maxGenotypeIndex));

    setVcfGenotypeFromGenotypeIndex(callerPloidy, maxGenotypeIndex, sampleInfo.maxGenotypeIndex);

    // set GQX
    sampleInfo.setGqx();

    // update homref prob for QUAL
    homRefLogProb += std::log(genotypePosterior[AG_GENOTYPE::HOMREF]);

    // set indelSampleInfo
    updateDiploidIndelLocusWithSampleInfo(opt, dopt, alleleGroup, sampleIndex, basecallBuffer, locus);

    // set haplotype ids and alt allele haplotype count ratio
    auto& maxGt = sampleInfo.max_gt();
    if (maxGt.getPloidy() >= 1)
    {
        auto allele0Index(maxGt.getAllele0Index());
        if (allele0Index > 0)
        {
            const auto& indelSampleData(alleleGroup.data(allele0Index-1).getSampleData(sampleIndex));
            maxGt.setAllele0HaplotypeId(indelSampleData.haplotypeId);
            maxGt.addAltAlleleHaplotypeCountRatio(indelSampleData.altAlleleHaplotypeCountRatio);
        }
    }
    if (maxGt.getPloidy() == 2)
    {
        auto allele1Index(maxGt.getAllele1Index());
        if (allele1Index > 0)
        {
            const auto& indelSampleData(alleleGroup.data(allele1Index-1).getSampleData(sampleIndex));
            maxGt.setAllele1HaplotypeId(indelSampleData.haplotypeId);
            maxGt.addAltAlleleHaplotypeCountRatio(indelSampleData.altAlleleHaplotypeCountRatio);
        }
    }
}



/// Determine if an allele group is reportable at \p pos
///
/// \return True if alleleGroup is not empty, and no allele positions are less than pos.
static
bool
isAlleleGroupReportable(
    const pos_t pos,
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup)
{
    const unsigned alleleGroupSize(alleleGroup.size());
    for (unsigned alleleIndex(0); alleleIndex < alleleGroupSize; ++alleleIndex)
    {
        const pos_t indelPos(alleleGroup.key(alleleIndex).pos);
        if (indelPos < pos)
        {
            return false;
        }
    }

    return (alleleGroupSize != 0);
}



/// Find all indel alleles at the current position which are either forced or viable candidates for
/// variant calling.
///
/// \param[out] orthogonalVariantAlleles All alleles at this position which meet acceptance criteria
static
void
getIndelAllelesAtPosition(
    const IndelBuffer& indelBuffer,
    const pos_t position,
    OrthogonalVariantAlleleCandidateGroup& orthogonalVariantAlleles)
{
    orthogonalVariantAlleles.clear();

    const unsigned sampleCount(indelBuffer.getSampleCount());

    auto it(indelBuffer.positionIterator(position));
    const auto it_end(indelBuffer.positionIterator(position + 1));
    for (; it != it_end; ++it)
    {
        const IndelKey& indelKey(it->first);

        if (indelKey.isMismatch()) continue;
        if (indelKey.is_breakpoint()) continue;

        const IndelData& indelData(getIndelData(it));

        const bool isForcedOutput(indelData.isForcedOutput);
        if (not isForcedOutput)
        {
            bool isZeroCoverage(true);
            for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
            {
                const IndelSampleData& indelSampleData(indelData.getSampleData(sampleIndex));

                if (not indelSampleData.read_path_lnp.empty())
                {
                    isZeroCoverage = false;
                    break;
                }
            }

            if (isZeroCoverage) continue;
            if (not indelBuffer.isCandidateIndel(indelKey, indelData)) continue;
        }

        // All alleles at the same position are automatically conflicting/orthogonal,
        // so they can be added to the orthogonal set without additional checks.
        orthogonalVariantAlleles.addVariantAllele(it);
    }
}



void
starling_pos_processor::
process_pos_indel_digt(const pos_t pos)
{
    // Check to see if we can call variant and forced alleles at this position.
    //
    // We may have already found this position's called and forced alleles while analyzing a locus positioned
    // downstream. If so, the entire locus is skipped over to prevent writing a conflicting interpretation of
    // the region covered by the downstream locus
    //
    if (pos <= _variantLocusAlreadyOutputToPos) return;


    // Reporting an extended indel locus is an all or nothing arrangement. An extended variant locus is composed
    // of a variant locus and all of the forced indels it overlaps. If any one member of this ensemble is
    // reportable, then the whole is treated as reportable. We rely on gvcf_writer to trim off any extra indels
    // left on the edges.
    bool isExtendedLocusReportable(is_pos_reportable(pos));


    // High-level summary of calling process:
    //
    // The locus calling process requires that we define groups of overlapping alleles to rank and then genotype.
    //
    // Overlapping alleles can be thought to form "conflict graphs", where an edge exists between two alleles
    // which cannot exist together on the same haplotype (called orthogonal alleles below). Without phasing
    // information, we can only (accurately) genotype among sets of alleles forming a clique in the graph.
    //
    // Given above constraint, we first identify all candidates alleles with a start position at the current
    // allele genotyper position (these form a clique by definition), and then greedily add the top-ranking
    // overlapping alleles with different start positions if they preserve the orthogonal clique relationship
    // of the set.
    //
    // Once we have the largest possible allele set, the reference is implicitly added and all alleles are
    // ranked. The top N are kept, N= ploidy. The reference is restored for the genotyping process if it is not
    // in the top N.
    //
    const unsigned sampleCount(getSampleCount());

    // First get all variant alleles at this position (including do-not-genotype alleles)
    OrthogonalVariantAlleleCandidateGroup positionVariantAlleles;
    getIndelAllelesAtPosition(getIndelBuffer(), pos, positionVariantAlleles);

    if (positionVariantAlleles.empty()) return;

    // Next remove the do-not-genotype alleles prior to computing all likelihoods
    // (these alleles will be added back in to the special forced allele pathway below)
    OrthogonalVariantAlleleCandidateGroup orthogonalVariantAlleles;
    {
        const unsigned positionAlleleGroupSize(positionVariantAlleles.size());
        for (unsigned alleleIndex(0); alleleIndex < positionAlleleGroupSize; ++alleleIndex)
        {
            const IndelData& indelData(positionVariantAlleles.data(alleleIndex));
            if (indelData.doNotGenotype) continue;

            orthogonalVariantAlleles.addVariantAllele(positionVariantAlleles.iter(alleleIndex));
        }
    }

#ifdef DEBUG_INDEL_OVERLAP
    log_os << "ZEBRA pos/pos-alleles: " << pos << " " << orthogonalVariantAlleles << "\n";
#endif

    // Determine ploidy for this locus in each sample
    //

    // Assume entire allele group is covered by one ploidy type per sample in nearly all cases,
    // in case of a conflict use the highest ploidy overlapped by the group.
    //
    const known_pos_range alleleGroupRange(positionVariantAlleles.getReferenceRange());

    std::vector<unsigned> groupLocusPloidy(sampleCount);
    for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
    {
        const unsigned groupLeftPloidy(get_ploidy(alleleGroupRange.begin_pos, sampleIndex));
        const unsigned groupRightPloidy(get_ploidy(alleleGroupRange.end_pos, sampleIndex));

        groupLocusPloidy[sampleIndex]=std::max(groupLeftPloidy,groupRightPloidy);
    }

    // A groupLocusPloidy value of 0 is treated as a special case, if this happens the
    // entire calling method reverts to a ploidy of 2 for the sample, as reflected in callerPloidy.
    //
    // The callerPloidy value will be used to generate the call information in this function,
    // and the groupLocusPloidy will be passed downstream to the gVCF writer.
    //
    std::vector<unsigned> callerPloidy(sampleCount);
    for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
    {
        callerPloidy[sampleIndex] = ((groupLocusPloidy[sampleIndex] == 0) ? 2 : groupLocusPloidy[sampleIndex]);
    }

    // This is set to true when at least one allele is reported based on variant evidence in at least one sample
    //
    // Even if this is not the case, other alleles may be reported because they are set to forcedOutput.
    bool isVariantLocusReported(false);

    // Track the top ranked variant alleles across all samples
    OrthogonalVariantAlleleCandidateGroup topVariantAlleleGroup;

    // The following subsection tracks data structures specific to variant allele calling phase of the process
    // which are not needed in the subsequent step where any remaining forcedOutput alleles are handled.
    if (! orthogonalVariantAlleles.empty())
    {
        // Track the top alt allele within each sample -- this is used as a temporary crutch to transfer the previous prior
        // calculation from single to multi-sample, and should be removed when a mature prior scheme is put in place
        std::vector<unsigned> topVariantAlleleIndexPerSample;

        // Rank input alleles to pick the top N per sample, where N=ploidy, and then aggregate and rank the top
        // per-sample alleles over all samples
        //
        // TODO: Note that the global allele ranking is unused by this particular function call, and could be optimized out.
        selectTopOrthogonalAllelesInAllSamples(
            sampleCount, callerPloidy, orthogonalVariantAlleles, topVariantAlleleGroup, topVariantAlleleIndexPerSample);

#ifdef DEBUG_INDEL_OVERLAP
        log_os << "ZEBRA pos/ranked-pos-alleles: " << pos << " " << topVariantAlleleGroup << "\n";
#endif

        // At this point topVariantAlleleGroup represents the alleles starting at the current position
        // which are best supported over all samples. Now we add conflicting alleles at other positions,
        // re-rank, and re-select the top alleles from the expanded set of candidates.
        if (not topVariantAlleleGroup.empty())
        {
            addAllelesAtOtherPositions(_ref, sampleCount, callerPloidy, pos,
                                       get_largest_total_indel_ref_span_per_read(),
                                       getIndelBuffer(), topVariantAlleleGroup, topVariantAlleleIndexPerSample);
        }

#ifdef DEBUG_INDEL_OVERLAP
        log_os << "ZEBRA pos/ranked-region-alleles: " << pos << " " << topVariantAlleleGroup << "\n";
#endif

        if (isAlleleGroupReportable(pos, topVariantAlleleGroup))
        {
            // The candidate variant alleles are eligible for reporting. At this point the alleles are
            // genotyped for each sample, the genotyped locus is tested to see if it meets emission criteria,
            // and if so, the locus information is pushed into the gVCF writer for reporting.
            //

            // parameter inputs if/when we wrap this as a function:
            static const bool isForcedOutput(false);
            static OrthogonalVariantAlleleCandidateGroup emptyGroup;

            // setup new indel locus:
            std::unique_ptr<GermlineIndelLocusInfo> locusPtr(
                new GermlineDiploidIndelLocusInfo(_dopt.gvcf, sampleCount));

            // cycle through variant alleles and add them to locus
            // (the locus interface requires that this is done before any other locus information is added):
            addIndelAllelesToLocus(topVariantAlleleGroup, isForcedOutput, *locusPtr);

            // Add sample-dependent info to locus:
            double homRefLogProb(0);
            ActiveRegionId activeRegionId(getActiveRegionDetector().getActiveRegionId(pos));
            for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
            {
                auto& sampleInfo(sample(sampleIndex));
                updateIndelLocusWithSampleInfo(
                    _opt, _dopt, topVariantAlleleGroup, topVariantAlleleIndexPerSample[sampleIndex], emptyGroup,
                    sampleInfo.sampleOptions,
                    callerPloidy[sampleIndex], groupLocusPloidy[sampleIndex], sampleIndex, activeRegionId, sampleInfo.basecallBuffer, *locusPtr,
                    homRefLogProb);
            }

            // Add sample-independent info to locus:
            locusPtr->anyVariantAlleleQuality = ln_error_prob_to_qphred(homRefLogProb);

            // get common prefix length
            if (locusPtr->getAltAlleleCount() > 1)
            {
                unsigned commonPrefixLength(getCommonPrefixLength(*locusPtr, _ref));
                if (commonPrefixLength > 0)
                {
                    locusPtr->setCommonPrefix(commonPrefixLength);
                }
            }

            if (locusPtr->isVariantLocus())
            {
                // The locus meets emission criteria, so send it into the gVCF writer for reporting.

                // Mark the end position of the range covered by this locus so that no interfering/duplicated
                // interpretations of this locus are written while analyzing downstream positions.
                //
                // The final position of the range goes one past the final locus position to represent interference of
                // adjacent indel alleles, for instance a 1D at position 10 should block any indel at position 11. This
                // one base offset does not need to be added because it's already built into the range objects
                // half-open interval format
                _variantLocusAlreadyOutputToPos = (locusPtr->range().end_pos());
                isVariantLocusReported = true;

                // STREL-392 check if the indel can be reported given this process' reporting restrictions
                //
                // Note two tricky details here:
                // 1. We need this check even though this function will not be called for "pos" values greater than
                //    the report range, because the new indel overlap resolution rules can lead to an indel locus
                //    begin_pos() which is higher than "pos", and thus potentially out of the report range.
                // 2. In order to better sync across independent processes, we need to call the enclosing function on a
                //    range of pos values before the report range, so that a value like "_variantLocusAlreadyOutputToPos"
                //    is in sync with the other process at the time we hit the first reportable position for this
                //    process, so we need to filter out indels found before report range begins
                //
                if (! isExtendedLocusReportable)
                {
                    isExtendedLocusReportable = is_pos_reportable(locusPtr->range().begin_pos());
                }

                if (isExtendedLocusReportable)
                {
                    _gvcfer->add_indel(std::move(locusPtr));
                }
            }
        }
    }

    //
    // Handle remaining forced variant output.
    //
    // In the section below the forced variants are handled using one of two approaches:
    //
    // 1) If a variant locus was reported at this position (ie. isVariantLocusReported is true)
    //
    // Handle all forced alleles from the current position to _variantLocusAlreadyOutputToPos
    //
    // 2) Otherwise...
    //
    // Handle all forced alleles at the current position only
    //

    // Determine which forced alleles have already been output via the standard variant calling process above.
    std::set<IndelKey> forcedAllelesAlreadyOutput;
    if (isVariantLocusReported)
    {
        // look through current topVariantAlleleGroup and note any forced output alleles already reported:
        const unsigned alleleGroupSize(topVariantAlleleGroup.size());
        for (unsigned alleleIndex(0); alleleIndex < alleleGroupSize; ++alleleIndex)
        {
            const IndelKey& indelKey(topVariantAlleleGroup.key(alleleIndex));
            const IndelData& indelData(topVariantAlleleGroup.data(alleleIndex));

            if (indelData.isForcedOutput)
            {
                forcedAllelesAlreadyOutput.insert(indelKey);
            }
        }
    }

    // score and report any remaining forced output alleles
    //
    {
        // Add all forced output alleles found at the current position and not already reported.
        //
        OrthogonalVariantAlleleCandidateGroup forcedOutputAlleleGroup;
        {
            const unsigned positionVariantAlleleCount(positionVariantAlleles.size());
            for (unsigned alleleIndex(0); alleleIndex < positionVariantAlleleCount; alleleIndex++)
            {
                const IndelKey& indelKey(positionVariantAlleles.key(alleleIndex));
                const IndelData& indelData(positionVariantAlleles.data(alleleIndex));

                const bool isForcedAndNotAlreadyOutput(indelData.isForcedOutput &&
                                                       (forcedAllelesAlreadyOutput.count(indelKey) == 0));

#ifdef FIXED_DONOTGENOTYPE_STATUS
                // All 'doNotGenotype' alleles should in theory have isForcedAndNotAlreadyOutput status, but this
                // assertion won't hold right now because complex indels can also come in directly from read
                // alignments. These will be marked as doNotGenotype even though they are not forced output
                //
                // TODO: fix doNotGenotype status such that this assertion can hold
                assert(isForcedAndNotAlreadyOutput || (! indelData.doNotGenotype));
#endif

                if (isForcedAndNotAlreadyOutput)
                {
                    forcedOutputAlleleGroup.addVariantAllele(positionVariantAlleles.iter(alleleIndex));
                }
            }
        }

        // If a variant locus has been reported at this position, then add forced output alleles
        // at all other positions covered by the reported locus
        if (isVariantLocusReported && (_variantLocusAlreadyOutputToPos > pos))
        {
            const auto& indelBuffer(getIndelBuffer());
            auto it(indelBuffer.positionIterator(pos+1));
            const auto it_end(indelBuffer.positionIterator(_variantLocusAlreadyOutputToPos+1));
            for (; it != it_end; ++it)
            {
                const IndelKey& indelKey(it->first);
                if (indelKey.is_breakpoint()) continue;
                if (indelKey.isMismatch()) continue;

                const IndelData& indelData(getIndelData(it));

                const bool isForcedAndNotAlreadyOutput(indelData.isForcedOutput &&
                                                       (forcedAllelesAlreadyOutput.count(indelKey) == 0));

#ifdef FIXED_DONOTGENOTYPE_STATUS
                // All 'doNotGenotype' alleles should in theory have isForcedAndNotAlreadyOutput status, see
                // comment in previous FIXED_DONOTGENOTYPE_STATUS block.
                assert(isForcedAndNotAlreadyOutput || (! indelData.doNotGenotype));
#endif

                if (isForcedAndNotAlreadyOutput)
                {
                    forcedOutputAlleleGroup.addVariantAllele(it);
                }
            }
        }

        static const bool isForcedOutput(true);

        // enumerate support for remaining forced output alleles compared to orthogonal genotyped variant alleles above
        const unsigned forcedOutputAlleleCount(forcedOutputAlleleGroup.size());
        for (unsigned forcedOutputAlleleIndex(0);
             forcedOutputAlleleIndex < forcedOutputAlleleCount; ++forcedOutputAlleleIndex)
        {
            // setup new indel locus:
            std::unique_ptr<GermlineIndelLocusInfo> locusPtr(new GermlineDiploidIndelLocusInfo(_dopt.gvcf, sampleCount));

            // fake an allele group with only the forced output allele so that we can output using
            // standard data structures
            OrthogonalVariantAlleleCandidateGroup fakeForcedOutputAlleleGroup;
            fakeForcedOutputAlleleGroup.addVariantAllele(forcedOutputAlleleGroup.alleles[forcedOutputAlleleIndex]);
            const unsigned fakeTopVariantAlleleIndexPerSample(0);

            // cycle through variant alleles and add them to locus (the locus interface requires that this is done first):
            addIndelAllelesToLocus(fakeForcedOutputAlleleGroup, isForcedOutput, *locusPtr);

            // add sample-dependent info:
            double homRefLogProb(0);
            ActiveRegionId activeRegionId(getActiveRegionDetector().getActiveRegionId(pos));
            for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
            {
                auto& sif(sample(sampleIndex));
                updateIndelLocusWithSampleInfo(
                    _opt, _dopt, fakeForcedOutputAlleleGroup, fakeTopVariantAlleleIndexPerSample, topVariantAlleleGroup,
                    sif.sampleOptions, callerPloidy[sampleIndex], groupLocusPloidy[sampleIndex], sampleIndex, activeRegionId,
                    sif.basecallBuffer, *locusPtr, homRefLogProb);
            }

            // add sample-independent info:
            locusPtr->anyVariantAlleleQuality = ln_error_prob_to_qphred(homRefLogProb);

            // STREL-392 check if the indel can be reported given this process' reporting restrictions
            //
            // see additional notes above for the previous add_indel call
            if (isExtendedLocusReportable)
            {
                // finished! send this locus down the pipe:
                _gvcfer->add_indel(std::move(locusPtr));
            }
        }
    }
}



/// \brief Translate the legacy indel allele support counts into the current AD/ADF/ADR reporting structure
///
/// The legacy structure was hard coded to only one indel allele. The new structure is designed to handle any number
/// of alts, as well as generalize over indel and SNV variant types. The new structure will be used to produce
/// AD/ADF/ADR counts in the VCF output and influences methods such as the LowDepth filter computation.
///
/// \param[in] legacyReportInfo - legacy indel allele support counts
/// \param[out] supportCounts - updated allele support counts
static
void
translateLegacyIndelSupportCounts(
    const AlleleSampleReportInfo& legacyReportInfo,
    LocusSupportingReadStats& supportCounts)
{
    // translating from legacy counting means there is only one alt allele by definition:
    static const unsigned nonRefAlleleCount(1);

    static const unsigned refAlleleIndex(0);
    static const unsigned altAlleleIndex(1);

    supportCounts.setAltCount(nonRefAlleleCount);

    supportCounts.fwdCounts.incrementAlleleCount(refAlleleIndex, legacyReportInfo.n_confident_ref_reads_fwd);
    supportCounts.fwdCounts.incrementAlleleCount(altAlleleIndex, legacyReportInfo.n_confident_indel_reads_fwd);
    supportCounts.fwdCounts.nonConfidentCount = legacyReportInfo.n_other_reads_fwd;

    supportCounts.revCounts.incrementAlleleCount(refAlleleIndex, legacyReportInfo.n_confident_ref_reads_rev);
    supportCounts.revCounts.incrementAlleleCount(altAlleleIndex, legacyReportInfo.n_confident_indel_reads_rev);
    supportCounts.revCounts.nonConfidentCount = legacyReportInfo.n_other_reads_rev;
}



/// Fill in all sample-specific indel locus info for the continuous-frequency calling case
///
static
void
updateContinuousIndelLocusWithSampleInfo(
    const starling_base_options& opt,
    const starling_deriv_options& dopt,
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    const unsigned sampleIndex,
    const pos_basecall_buffer& basecallBuffer,
    GermlineIndelLocusInfo& locus)
{
    auto& sampleInfo(locus.getSample(sampleIndex));
    sampleInfo.setPloidy(-1);

    assert(alleleGroup.size() == 1);

    static const bool is_use_alt_indel(true);
    updateIndelSampleInfo(opt,dopt, alleleGroup, sampleIndex, basecallBuffer, is_use_alt_indel, locus);

    const GermlineIndelSampleInfo& indelSampleInfo(locus.getIndelSample(sampleIndex));

    // translate the legacy read support counts into the current AD/ADF/ADR reporting structure
    translateLegacyIndelSupportCounts(indelSampleInfo.legacyReportInfo, sampleInfo.supportCounts);

    sampleInfo.gqx = sampleInfo.genotypeQualityPolymorphic =
                         starling_continuous_variant_caller::getAlleleSequencingErrorQscore(
                             indelSampleInfo.legacyReportInfo.n_confident_indel_reads,
                             indelSampleInfo.legacyReportInfo.total_confident_reads(),
                             opt.continuousFrequencyCallerExpectedObservationQuality,
                             opt.continuousFrequencyCallerMaxQscore);

    // use diploid gt codes as a convenient way to summarize the continuous variant calls:
    static const VcfGenotype hetGtIndex(0,1);
    static const VcfGenotype homGtIndex(1,1);

    const bool isHetLike(indelSampleInfo.alleleFrequency() < (1 - opt.min_het_vf));

    sampleInfo.maxGenotypeIndexPolymorphic = (isHetLike ? hetGtIndex : homGtIndex);
}



void
starling_pos_processor::
process_pos_indel_continuous(const pos_t pos)
{
    const unsigned sampleCount(getSampleCount());

    OrthogonalVariantAlleleCandidateGroup orthogonalVariantAlleles;
    getIndelAllelesAtPosition(getIndelBuffer(), pos, orthogonalVariantAlleles);

    // In continuous calling mode, each non-reference allele is handled independently:
    //
    const unsigned nonrefAlleleCount(orthogonalVariantAlleles.size());
    for (unsigned nonrefAlleleIndex(0); nonrefAlleleIndex<nonrefAlleleCount; nonrefAlleleIndex++)
    {
        const auto nonrefAlleleIter(orthogonalVariantAlleles.iter(nonrefAlleleIndex));
        const bool isForcedOutput(orthogonalVariantAlleles.data(nonrefAlleleIndex).isForcedOutput);

        OrthogonalVariantAlleleCandidateGroup topVariantAlleleGroup;
        topVariantAlleleGroup.addVariantAllele(nonrefAlleleIter);

        // setup new indel locus:
        std::unique_ptr<GermlineContinuousIndelLocusInfo> locusPtr(new GermlineContinuousIndelLocusInfo(sampleCount));

        // cycle through variant alleles and add them to locus (the locus interface requires that this is done first):
        addIndelAllelesToLocus(topVariantAlleleGroup, isForcedOutput, *locusPtr);

        // add sample-dependent info:
        for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
        {
            updateContinuousIndelLocusWithSampleInfo(
                _opt, _dopt, topVariantAlleleGroup, sampleIndex, sample(sampleIndex).basecallBuffer, *locusPtr);
        }

        // add sample-independent info:
        {
            int anyVariantAlleleQuality = 0;
            for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
            {
                auto& sampleInfo(locusPtr->getSample(sampleIndex));
                anyVariantAlleleQuality = std::max(anyVariantAlleleQuality,
                                                   sampleInfo.genotypeQualityPolymorphic);
            }
            locusPtr->anyVariantAlleleQuality = anyVariantAlleleQuality;
        }

        // determine if this locus (in continuous case locus TEMPORARILY means "allele") is printable:
        bool isReportableLocus(isForcedOutput);
        for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
        {
            const auto& indelSampleInfo(locusPtr->getIndelSample(sampleIndex));
            if (not isReportableLocus)
            {
                const double alleleFrequency(indelSampleInfo.alleleFrequency());
                if ((indelSampleInfo.legacyReportInfo.n_confident_indel_reads > 0) and
                    (alleleFrequency > _opt.min_het_vf))
                {
                    isReportableLocus = true;
                }
            }
        }
        if (not isReportableLocus) continue;

        _gvcfer->add_indel(std::move(locusPtr));
    }
}

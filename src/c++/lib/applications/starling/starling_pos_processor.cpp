//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2017 Illumina, Inc.
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
    : base_t(opt, dopt, ref, fileStreams, opt.alignFileOpt.alignmentFilename.size(), statsManager),
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
    _forcedAllelesAlreadyOutput.clear();
    _active_region_detector->clear();
}



void
starling_pos_processor::
process_pos_snp(const pos_t pos)
{
    try
    {
        const unsigned sampleCount(getSampleCount());

        const bool isForcedOutput(is_forced_output_pos(pos));

        // the second term in is_skippable below forces sites to go through the pipeline
        // if phaser has put a hold on buffer cleanup. This ensures that the phaser will be turned back off
        //
        // TODO: there must be a way to force correct usage into the phaser's API instead of requiring this brittle hack
        const bool isSkippable(!isForcedOutput);

        if (isSkippable)
        {
            bool isZeroCoverage(true);
            for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
            {
                const sample_info& sif(sample(sampleIndex));
                const CleanedPileup& cpi(sif.cpi);
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
            _pileupCleaner.CleanPileupErrorProb(sample(sampleIndex).cpi);
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



/// setup siteSampleInfo assuming that corresponding sampleInfo has already been initialized
static
void
updateSiteSampleInfo(
    const starling_base_options& opt,
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

    siteSampleInfo.n_used_calls=cpi.n_used_calls();
    siteSampleInfo.n_unused_calls=cpi.n_unused_calls();

    // MQ is computed/reported whether EVS features are needed or not, it is also used by EVS
    const auto& pi(cpi.rawPileup());
    siteSampleInfo.mapqTracker = pi.mapqTracker;

    siteSampleInfo.strandBias = strandBias;

    /// add EVS feature info
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
    const extended_pos_info& good_epi(sif.cpi.getExtendedPosInfo());
    dgt.ploidy=ploidy;
    dopt.pdcaller().position_snp_call_pprob_digt(
        opt, good_epi, dgt, opt.is_all_sites());
}



/// translate a base index into the an allele index for the locus
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
        oss << "ERROR: Can't find allele " << id_to_base(baseIndex) << " expected to be present in site locus\n";
        BOOST_THROW_EXCEPTION(LogicException(oss.str()));
    }

    return alleleIndex;
}



/// translate older starling "4-allele" genotype index to current VcfGenotype associated
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
    const starling_base_options& opt,
    const starling_pos_processor::sample_info& sif,
    const unsigned callerPloidy,
    const unsigned groupLocusPloidy,
    const diploid_genotype& dgt,
    const unsigned sampleIndex,
    const CandidateSnvBuffer& candidateSnvBuffer,
    GermlineDiploidSiteLocusInfo& locus,
    double& homRefLogProb)
{
    auto& sampleInfo(locus.getSample(sampleIndex));
    sampleInfo.setPloidy(callerPloidy);

    bool isOverlappingHomAltDeletion(false);
    if (groupLocusPloidy == 0)
    {
        const int locusPloidyAdjustment(sif.cpi.rawPileup().spanningIndelPloidyModification);
        isOverlappingHomAltDeletion=(locusPloidyAdjustment < 0);
    }

    const CleanedPileup& cpi(sif.cpi);

    if (cpi.n_used_calls() != 0)
    {
        // the principle of this filter is that there's supposed to be no coverage here
        // we make an exception for sites inside of homalt deletions, maybe we shouldn't?
        if ((groupLocusPloidy == 0) and (not isOverlappingHomAltDeletion))
        {
            sampleInfo.setPloidyConflict();
        }
    }

    if     (locus.isRefUnknown() or (cpi.n_used_calls() == 0) or isOverlappingHomAltDeletion)
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
                if (call.get_qscore() < opt.used_allele_count_min_qscore) continue;
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
    if (maxGt.getPloidy() == 2)
    {
        // diploid
        auto allele0Index(maxGt.getAllele0Index());
        auto allele1Index(maxGt.getAllele1Index());

        if (allele0Index > 0)
            maxGt.setAllele0HaplotypeId(
                candidateSnvBuffer.getHaplotypeId(
                    sampleIndex, locus.pos, locus.getSiteAlleles()[allele0Index-1].baseIndex)
            );
        if (allele1Index > 0)
            maxGt.setAllele1HaplotypeId(
                candidateSnvBuffer.getHaplotypeId(
                    sampleIndex, locus.pos, locus.getSiteAlleles()[allele1Index-1].baseIndex)
            );
    }
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
        const auto& cpi(sample(sampleIndex).cpi);
        const auto& good_pi(cpi.cleanedPileup());
        good_pi.get_known_counts(sampleBaseCounts, _opt.used_allele_count_min_qscore);

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
        const int locusPloidyAdjustment(sample(sampleIndex).cpi.rawPileup().spanningIndelPloidyModification);
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
    auto activeRegionId(_active_region_detector->getActiveRegionId(pos));
    std::unique_ptr<GermlineDiploidSiteLocusInfo> locusPtr(new GermlineDiploidSiteLocusInfo(_dopt.gvcf, sampleCount, activeRegionId, pos, refBaseIndex, isForcedOutput));


    // add all candidate alternate alleles:
    for (const auto baseId : altAlleles)
    {
        locusPtr->addAltSiteAllele(static_cast<BASE_ID::index_t>(baseId));
    }

    double homRefLogProb(0);
    for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
    {
        updateSnvLocusWithSampleInfo(
            _opt, sample(sampleIndex), callerPloidy[sampleIndex], groupLocusPloidy[sampleIndex],
            allDgt[sampleIndex], sampleIndex, _active_region_detector->getCandidateSnvBuffer(), *locusPtr, homRefLogProb);
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
    const starling_base_options& opt,
    starling_pos_processor::sample_info& sif,
    const unsigned sampleIndex,
    std::vector<StrandBiasCounts>& strandBiasCounts,
    GermlineContinuousSiteLocusInfo& locus)
{
    const auto& siteAlleles(locus.getSiteAlleles());
    const bool isRefAllele(siteAlleles.empty());

    const uint8_t altAlleleCount(siteAlleles.size());
    const uint8_t fullAlleleCount(altAlleleCount+1);

    const CleanedPileup& cpi(sif.cpi);
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
            if (call.get_qscore()<opt.used_allele_count_min_qscore) continue;
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
                             starling_continuous_variant_caller::poisson_qscore(
                                 continuousSiteSampleInfo.continuousAlleleDepth,
                                 continuousSiteSampleInfo.continuousTotalDepth,
                                 (unsigned) opt.continuousSiteCallerAverageQuality, 40);
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
        // regions of the genome, we need to go through a dummy calling operation for approx 'max_indel_size' bases
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
static unsigned getCommonPrefixLength(
    const GermlineIndelLocusInfo& locus,
    const reference_contig_segment& ref)
{
    const auto& locusRange(locus.range());

    unsigned minLenCommonPrefix(locusRange.size());
    for (auto& indelAllele : locus.getIndelAlleles())
    {
        auto& indelKey = indelAllele.indelKey;

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
                    break;
                ++lenCommonPrefix;
            }
        }
        if (lenCommonPrefix < minLenCommonPrefix)
            minLenCommonPrefix = lenCommonPrefix;
    }

    return minLenCommonPrefix;
}


/// setup indelSampleInfo assuming that corresponding sampleInfo has already been initialized
///
/// \param basecallBuffer[inout] read previous base depth and mqpq values from basecallbuffer, and also
///                              set ploidy modifications on any bases which the indel overlaps
static
void
updateIndelSampleInfo(
    const starling_base_options& opt,
    const starling_deriv_options& dopt,
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    const unsigned sampleIndex,
    pos_basecall_buffer& basecallBuffer,
    GermlineIndelLocusInfo& locus)
{
    GermlineIndelSampleInfo indelSampleInfo;

    const LocusSampleInfo& sampleInfo(locus.getSample(sampleIndex));
    const unsigned callerPloidy(sampleInfo.getPloidy().getPloidy());

    // set site ploidy modifications:
    const auto& range(locus.range());
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

    {
        // get various indel stats from the pileup:
        pos_t pileupPos(range.begin_pos()-1);
        const IndelKey& indelKey0(alleleGroup.key(0));
        if (indelKey0.type == INDEL::BP_RIGHT) pileupPos=range.end_pos();
        const snp_pos_info& spi(basecallBuffer.get_pos(pileupPos));
        indelSampleInfo.tier1Depth=spi.calls.size();
        indelSampleInfo.mapqTracker=spi.mapqTracker;
    }

    /// TODO STREL-125 deprecated
    // add misc sample info from legacy sample indel report:
    {
        static const bool is_tier2_pass(false);
        static const bool is_use_alt_indel(false);

        /// TODO STREL-125 legacy structure assumes single indel allele, get rid of this....
        const IndelKey& indelKey(alleleGroup.key(0));
        const IndelData& indelData(alleleGroup.data(0));
        const IndelSampleData& indelSampleData(indelData.getSampleData(sampleIndex));
        getAlleleSampleReportInfo(opt, dopt, indelKey, indelSampleData, basecallBuffer,
                                  is_tier2_pass, is_use_alt_indel, indelSampleInfo.legacyReportInfo);
    }

    locus.setIndelSampleInfo(sampleIndex, indelSampleInfo);
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
/// \param basecallBuffer[inout] see updateIndelSampleInfo
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
    const unsigned patternRepeatCount=std::max(1u,indelReportInfo.ref_repeat_count);

    const ContextGenotypePriors& genotypePriors(dopt.getIndelGenotypePriors().getContextSpecificPriorSet(patternRepeatCount));

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
    normalize_ln_distro(std::begin(genotypePosterior), std::end(genotypePosterior), maxGenotypeIndex);

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

    normalize_ln_distro(std::begin(genotypePosterior), std::end(genotypePosterior), maxGenotypeIndex);

    sampleInfo.genotypeQuality = error_prob_to_qphred(prob_comp(std::begin(genotypePosterior), std::end(genotypePosterior), maxGenotypeIndex));

    setVcfGenotypeFromGenotypeIndex(callerPloidy, maxGenotypeIndex, sampleInfo.maxGenotypeIndex);

    // set GQX
    sampleInfo.setGqx();

    // update homref prob for QUAL
    homRefLogProb += std::log(genotypePosterior[AG_GENOTYPE::HOMREF]);

    // set indelSampleInfo
    updateIndelSampleInfo(opt, dopt, alleleGroup, sampleIndex, basecallBuffer, locus);

    // set haplotype ids
    auto& maxGt = sampleInfo.max_gt();
    if (maxGt.getPloidy() == 2)
    {
        // diploid
        auto allele0Index(maxGt.getAllele0Index());
        auto allele1Index(maxGt.getAllele1Index());

        if (allele0Index > 0)
            maxGt.setAllele0HaplotypeId(alleleGroup.data(allele0Index-1).getSampleData(sampleIndex).haplotypeId);
        if (allele1Index > 0)
            maxGt.setAllele1HaplotypeId(alleleGroup.data(allele1Index-1).getSampleData(sampleIndex).haplotypeId);
    }
}



/// determine if an allele group is reportable
///
/// \return true if alleles exist and no alleles are less than pos.
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



void
starling_pos_processor::
process_pos_indel_digt(const pos_t pos)
{
    const unsigned sampleCount(getSampleCount());

    auto it(getIndelBuffer().positionIterator(pos));
    const auto it_end(getIndelBuffer().positionIterator(pos + 1));

    // define groups of overlapping alleles to rank and then genotype.
    //
    // overlapping alleles can be thought to form "conflict graphs", where an edge exists between two alleles
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
    OrthogonalVariantAlleleCandidateGroup orthogonalVariantAlleles;
    for (; it!=it_end; ++it)
    {
        const IndelKey& indelKey(it->first);
        const IndelData& indelData(getIndelData(it));

        if (indelKey.is_breakpoint()) continue;

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
            if (not getIndelBuffer().isCandidateIndel(indelKey, indelData)) continue;
        }

        // all alleles at the same position are automatically conflicting/orthogonal:
        orthogonalVariantAlleles.addVariantAllele(it);
    }

    if (orthogonalVariantAlleles.empty()) return;

#ifdef DEBUG_INDEL_OVERLAP
    log_os << "ZEBRA pos/pos-alleles: " << pos << " " << orthogonalVariantAlleles << "\n";
#endif

    // determine ploidy for this locus in each sample
    //
    std::vector<unsigned> groupLocusPloidy(sampleCount);
    for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
    {
        // Assume entire allele group is covered by one ploidy type per sample in nearly all cases,
        // in case of a conflict use the highest ploidy overlapped by the group.
        //
        const known_pos_range alleleGroupRange(orthogonalVariantAlleles.getReferenceRange());
        const unsigned groupLeftPloidy(get_ploidy(alleleGroupRange.begin_pos, sampleIndex));
        const unsigned groupRightPloidy(get_ploidy(alleleGroupRange.end_pos, sampleIndex));

        groupLocusPloidy[sampleIndex]=std::max(groupLeftPloidy,groupRightPloidy);
    }

    // groupLocusPloidy of 0 is treated as a special case, if this happens the
    // entire calling method reverts to a ploidy of 2 for the sample, but the
    // locus ploidy is passed into the gVCF writer as 0. The gVCF writer can
    // decide what to do with this information from there.
    //
    std::vector<unsigned> callerPloidy(sampleCount);
    for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
    {
        callerPloidy[sampleIndex] = ((groupLocusPloidy[sampleIndex] == 0) ? 2 : groupLocusPloidy[sampleIndex]);
    }

    bool isReportedLocus(false);
    OrthogonalVariantAlleleCandidateGroup topVariantAlleleGroup;

    auto activeRegionId(_active_region_detector->getActiveRegionId(pos));

    // now check to see if we can call variant alleles at this position, we may have already found this positions alleles while
    // analyzing a locus positioned downstream. If so, skip ahead to handle the forced output alleles only:
    if (_variantLocusAlreadyOutputToPos <= pos)
    {
        // track top alt allele within each sample -- this is used as a temporary crutch to transfer the previous prior
        // calculation from single to multi-sample, and should be removed when a mature prior scheme is put in place
        std::vector<unsigned> topVariantAlleleIndexPerSample(sampleCount);

        // rank input alleles to pick the top N, N=ploidy, per sample, and aggregate/rank these
        // over all samples
        selectTopOrthogonalAllelesInAllSamples(
            sampleCount, callerPloidy, orthogonalVariantAlleles, topVariantAlleleGroup, topVariantAlleleIndexPerSample);

#ifdef DEBUG_INDEL_OVERLAP
        log_os << "ZEBRA pos/ranked-pos-alleles: " << pos << " " << topVariantAlleleGroup << "\n";
#endif

        // At this point topVariantAlleleGroup represents the best alleles which
        // start at the current position (over all samples). Now we add conflicting
        // alleles at other positions and re-rank, re-select the top alleles again:
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
            // genotype and report topVariantAlleleGroup
            //

            // parameter inputs if/when we wrap this as a function:
            static const bool isForcedOutput(false);
            static OrthogonalVariantAlleleCandidateGroup emptyGroup;

            // setup new indel locus:
            std::unique_ptr<GermlineIndelLocusInfo> locusPtr(
                new GermlineDiploidIndelLocusInfo(_dopt.gvcf, sampleCount, activeRegionId));

            // cycle through variant alleles and add them to locus (the locus interface requires that this is done first):
            addIndelAllelesToLocus(topVariantAlleleGroup, isForcedOutput, *locusPtr);

            // add sample-dependent info:
            double homRefLogProb(0);
            for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
            {
                auto& sif(sample(sampleIndex));
                updateIndelLocusWithSampleInfo(
                    _opt, _dopt, topVariantAlleleGroup, topVariantAlleleIndexPerSample[sampleIndex], emptyGroup,
                    sif.sample_opt,
                    callerPloidy[sampleIndex], groupLocusPloidy[sampleIndex], sampleIndex, sif.bc_buff, *locusPtr,
                    homRefLogProb);
            }

            // add sample-independent info:
            locusPtr->anyVariantAlleleQuality = ln_error_prob_to_qphred(homRefLogProb);

            // get common prefix length

            if (locusPtr->getAltAlleleCount() > 1)
            {
                unsigned commonPrefixLength(getCommonPrefixLength(*locusPtr, _ref));
                if (commonPrefixLength > 0)
                    locusPtr->setCommonPrefix(commonPrefixLength);
            }

            if (isForcedOutput or locusPtr->isVariantLocus())
            {
                // finished! send this locus down the pipe:

                // expand the end range of the locus by one to represent adjacent indel interference, for instance a
                // 1D at position 10 should block any indel at position 11
                _variantLocusAlreadyOutputToPos = (locusPtr->range().end_pos() + 1);
                isReportedLocus = true;

                // STREL-392 check if the indel can be reported given this process' reporting restrictions
                //
                // Note two tricky details here:
                // 1. we need this check even though this function will not be called for "pos" values greater than
                //    the report range, because the new indel overlap resolution rules can lead to an indel locus
                //    begin_pos() which is higher than "pos", and thus potentially out of the report range.
                // 2. In order to better sync across independent processes, we need to call the enclosing function on a
                //    range of pos values before the report range, so that values like "_variantLocusAlreadyOutputToPos"
                //    and "_forcedAllelesAlreadyOutput" are in sync with the other process at the time we hit the first
                //    reportable position for this process, so we need to filter out indels found before report range
                //    begins
                //
                if (is_pos_reportable(locusPtr->range().begin_pos()))
                {
                    _gvcfer->add_indel(std::move(locusPtr));
                }
            }
        }
    }

    // update data structure to track which forced alleles have already been output ahead of the current position
    {
        // we never track anything below the current position
        const auto endIter(_forcedAllelesAlreadyOutput.lower_bound(IndelKey(pos)));
        _forcedAllelesAlreadyOutput.erase(_forcedAllelesAlreadyOutput.begin(),endIter);

        if (isReportedLocus)
        {
            // look through current topVariantAlleleGroup and note any forced output alleles already reported:
            const unsigned alleleGroupSize(topVariantAlleleGroup.size());
            for (unsigned alleleIndex(0); alleleIndex < alleleGroupSize; ++alleleIndex)
            {
                const IndelKey& indelKey(topVariantAlleleGroup.key(alleleIndex));
                const IndelData& indelData(topVariantAlleleGroup.data(alleleIndex));

                if (indelData.isForcedOutput)
                {
                    _forcedAllelesAlreadyOutput.insert(indelKey);
                }
            }
        }
    }

    // score and report any remaining forced output alleles
    //
    {
        // track all forced output alleles in a separate group (even if they go into topVariant group)
        // to ensure that these are output even if not included in the most likely genotype for any sample
        //
        OrthogonalVariantAlleleCandidateGroup forcedOutputAlleleGroup;
        {
            const unsigned orthogonalVariantAlleleCount(orthogonalVariantAlleles.size());
            for (unsigned alleleIndex(0); alleleIndex < orthogonalVariantAlleleCount; alleleIndex++)
            {
                const IndelKey& indelKey(orthogonalVariantAlleles.key(alleleIndex));
                const IndelData& indelData(orthogonalVariantAlleles.data(alleleIndex));
                if (indelData.isForcedOutput and (_forcedAllelesAlreadyOutput.count(indelKey) == 0))
                {
                    forcedOutputAlleleGroup.addVariantAllele(orthogonalVariantAlleles.iter(alleleIndex));
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
            std::unique_ptr<GermlineIndelLocusInfo> locusPtr(new GermlineDiploidIndelLocusInfo(_dopt.gvcf, sampleCount, activeRegionId));

            // fake an allele group with only the forced output allele so that we can output using
            // standard data structures
            OrthogonalVariantAlleleCandidateGroup fakeForcedOutputAlleleGroup;
            fakeForcedOutputAlleleGroup.addVariantAllele(forcedOutputAlleleGroup.alleles[forcedOutputAlleleIndex]);
            const unsigned fakeTopVariantAlleleIndexPerSample(0);

            // cycle through variant alleles and add them to locus (the locus interface requires that this is done first):
            addIndelAllelesToLocus(fakeForcedOutputAlleleGroup, isForcedOutput, *locusPtr);

            // add sample-dependent info:
            double homRefLogProb(0);
            for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
            {
                auto& sif(sample(sampleIndex));
                updateIndelLocusWithSampleInfo(
                    _opt, _dopt, fakeForcedOutputAlleleGroup, fakeTopVariantAlleleIndexPerSample, topVariantAlleleGroup,
                    sif.sample_opt, callerPloidy[sampleIndex], groupLocusPloidy[sampleIndex], sampleIndex,
                    sif.bc_buff, *locusPtr, homRefLogProb);
            }

            // add sample-independent info:
            locusPtr->anyVariantAlleleQuality = ln_error_prob_to_qphred(homRefLogProb);

            // STREL-392 check if the indel can be reported given this process' reporting restrictions
            //
            // see additional notes above for the previous add_indel call
            if (is_pos_reportable(locusPtr->range().begin_pos()))
            {
                // finished! send this locus down the pipe:
                _gvcfer->add_indel(std::move(locusPtr));
            }
        }
    }
}



/// fill in all sample-specific locus info
///
static
void
updateContinuousIndelLocusWithSampleInfo(
    const starling_base_options& opt,
    const starling_deriv_options& dopt,
    const unsigned sampleIndex,
    const pos_basecall_buffer& bc_buff,
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    GermlineIndelLocusInfo& locus)
{
    static const bool is_tier2_pass(false);
    static const bool is_use_alt_indel(true);

    auto& sampleInfo(locus.getSample(sampleIndex));
    sampleInfo.setPloidy(-1);

    // tmp until new count structures are in place:
    assert(alleleGroup.size() == 1);
    const auto& indelKey(alleleGroup.key(0));
    const auto& indelData(alleleGroup.data(0));

    const IndelSampleData& indelSampleData(indelData.getSampleData(sampleIndex));

    GermlineIndelSampleInfo indelSampleInfo;
    getAlleleSampleReportInfo(opt, dopt, indelKey, indelSampleData, bc_buff, is_tier2_pass,
                              is_use_alt_indel, indelSampleInfo.legacyReportInfo);
    locus.setIndelSampleInfo(sampleIndex, indelSampleInfo);

    sampleInfo.gqx = sampleInfo.genotypeQualityPolymorphic =
                         starling_continuous_variant_caller::poisson_qscore(
                             indelSampleInfo.legacyReportInfo.n_confident_indel_reads,
                             indelSampleInfo.legacyReportInfo.total_confident_reads(),
                             (unsigned) opt.continuousSiteCallerAverageQuality, 40);

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

    auto it(getIndelBuffer().positionIterator(pos));
    const auto it_end(getIndelBuffer().positionIterator(pos + 1));
    for (; it!=it_end; ++it)
    {
        const IndelKey& indelKey(it->first);
        const IndelData& indelData(getIndelData(it));

        if (indelKey.is_breakpoint()) continue;

        const bool isForcedOutput(indelData.isForcedOutput);

        if (! isForcedOutput)
        {
            bool isZeroCoverage(true);
            for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
            {
                const IndelSampleData& indelSampleData(indelData.getSampleData(sampleIndex));
                if (not indelSampleData.read_path_lnp.empty())
                {
                    isZeroCoverage = false;
                    break;
                }
            }

            if (isZeroCoverage) continue;
            if (!getIndelBuffer().isCandidateIndel(indelKey, indelData)) continue;
        }


        // the way things are handled in continuous mode right now we only add one allele per locus
        OrthogonalVariantAlleleCandidateGroup topVariantAlleleGroup;
        topVariantAlleleGroup.alleles.push_back(it);

        // setup new indel locus:
        std::unique_ptr<GermlineContinuousIndelLocusInfo> locusPtr(new GermlineContinuousIndelLocusInfo(sampleCount));

        // cycle through variant alleles and add them to locus (the locus interface requires that this is done first):
        addIndelAllelesToLocus(topVariantAlleleGroup, isForcedOutput, *locusPtr);

        // add sample-dependent info:
        for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
        {
            updateContinuousIndelLocusWithSampleInfo(
                _opt, _dopt, sampleIndex, sample(sampleIndex).bc_buff, topVariantAlleleGroup, *locusPtr);
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

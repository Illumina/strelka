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

#include "starling_pos_processor.hh"

#include "starling_continuous_variant_caller.hh"
#include "blt_common/ref_context.hh"
#include "blt_util/log.hh"
#include "blt_util/prob_util.hh"
#include "blt_util/sort_util.hh"
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
    const starling_streams& streams)
    : base_t(opt,dopt,ref,streams, opt.alignFileOpt.alignmentFilename.size()),
      _opt(opt),
      _dopt(dopt),
      _streams(streams)
{
    static const unsigned sampleId(0);

    assert(_streams.getSampleNames().size() == getSampleCount());

    // setup gvcf aggregator
    if (_opt.gvcf.is_gvcf_output())
    {
        _gvcfer.reset(new gvcf_aggregator(
                          _opt,_dopt,ref,_nocompress_regions,
                          _streams.getSampleNames(), _streams.gvcf_osptr(),
                          sample(sampleId).bc_buff));
    }

    // setup indel buffer:
    {
        double maxIndelCandidateDepthSumOverNormalSamples(-1.);

        if (dopt.gvcf.is_max_depth())
        {
            if (opt.max_candidate_indel_depth_factor > 0.)
            {
                maxIndelCandidateDepthSumOverNormalSamples = (opt.max_candidate_indel_depth_factor * dopt.gvcf.max_depth);
            }
        }

        if (opt.max_candidate_indel_depth > 0.)
        {
            if (maxIndelCandidateDepthSumOverNormalSamples > 0.)
            {
                maxIndelCandidateDepthSumOverNormalSamples = std::min(maxIndelCandidateDepthSumOverNormalSamples,static_cast<double>(opt.max_candidate_indel_depth));
            }
            else
            {
                maxIndelCandidateDepthSumOverNormalSamples = opt.max_candidate_indel_depth;
            }
        }

        getIndelBuffer().setMaxCandidateDepth(maxIndelCandidateDepthSumOverNormalSamples);

        const unsigned sampleCount(getSampleCount());
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
insert_nocompress_region(
    const known_pos_range2& range)
{
    _stageman.validate_new_pos_value(range.begin_pos(),STAGE::READ_BUFFER);
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
        const bool isSkippable(! (isForcedOutput || is_save_pileup_buffer()));

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
    GermlineSiteLocusInfo& locus)
{
    const snp_pos_info& good_pi(cpi.cleanedPileup());

    GermlineSiteSampleInfo siteSampleInfo;

    siteSampleInfo.spanningDeletionReadCount = good_pi.spanningDeletionReadCount;

    siteSampleInfo.n_used_calls=cpi.n_used_calls();
    siteSampleInfo.n_unused_calls=cpi.n_unused_calls();

    good_pi.get_known_counts(siteSampleInfo.fwd_counts, opt.used_allele_count_min_qscore, true);
    good_pi.get_known_counts(siteSampleInfo.rev_counts, opt.used_allele_count_min_qscore, false);

    // MQ is computed/reported whether EVS features are needed or not, it is also used by EVS
    const auto& pi(cpi.rawPileup());
    siteSampleInfo.mapqTracker = pi.mapqTracker;

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
        if (allele.baseIndex == baseIndex)
        {
            return alleleIndex;
        }
        alleleIndex++;
    }

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
    GermlineDiploidSiteLocusInfo& locus,
    double& homRefLogProb)
{
    auto& sampleInfo(locus.getSample(sampleIndex));
    sampleInfo.setPloidy(callerPloidy);
    if (callerPloidy != groupLocusPloidy)
    {
        sampleInfo.setPloidyConflict();
    }

    const CleanedPileup& cpi(sif.cpi);
    //const extended_pos_info& good_epi(cpi.getExtendedPosInfo());

    sampleInfo.genotypeQuality = dgt.genome.max_gt_qphred;
    translateDigtToVcfGenotype(callerPloidy, dgt.genome.max_gt, locus, sampleInfo.maxGenotypeIndex);

    if     (locus.isRefUnknown() || (cpi.n_used_calls()==0))
    {
        sampleInfo.genotypeQualityPolymorphic=0;
        sampleInfo.maxGenotypeIndexPolymorphic.setGenotypeFromAlleleIndices();
        sampleInfo.gqx=0;
    }
    else
    {
        sampleInfo.genotypeQualityPolymorphic = dgt.poly.max_gt_qphred;
        translateDigtToVcfGenotype(callerPloidy, dgt.poly.max_gt, locus, sampleInfo.maxGenotypeIndexPolymorphic);

        sampleInfo.setGqx();
    }

    // set PL values:
    const auto& siteAlleles(locus.getSiteAlleles());
    const uint8_t nonRefAlleleCount(siteAlleles.size());
    const bool isAltAlleles(nonRefAlleleCount>0);
    if (isAltAlleles)
    {
        const uint8_t fullAlleleCount(nonRefAlleleCount+1);

        // number of PL fields required:
        //const unsigned genotypeCount(VcfGenotypeUtil::getGenotypeCount(callerPloidy, fullAlleleCount));
        sampleInfo.genotypePhredLoghood.setPloidy(callerPloidy);

        if (callerPloidy == 1)
        {
            for (unsigned allele0Index(0); allele0Index<=fullAlleleCount; ++allele0Index)
            {
                const uint8_t base0Index(siteAlleles[allele0Index].baseIndex);
                DIGT::get_gt_with_alleles(base0Index,base0Index);
                sampleInfo.genotypePhredLoghood.getGenotypeLikelihood(allele0Index) =
                    dgt.phredLoghood[base0Index];
            }
        }
        else if(callerPloidy == 2)
        {
            for (unsigned allele1Index(0); allele1Index < fullAlleleCount; ++allele1Index)
            {
                for (unsigned allele0Index(0); allele0Index <= allele1Index; ++allele0Index)
                {
                    const uint8_t base0Index(siteAlleles[allele0Index].baseIndex);
                    const uint8_t base1Index(siteAlleles[allele1Index].baseIndex);
                    const unsigned digtGenotypeIndex(DIGT::get_gt_with_alleles(base0Index,base1Index));
                    sampleInfo.genotypePhredLoghood.getGenotypeLikelihood(allele0Index,allele1Index) =
                        dgt.phredLoghood[digtGenotypeIndex];
                }
            }
        }
        else
        {
            assert(false and "Unexpected ploidy");
        }
    }

    // update homref prob for QUAL
    homRefLogProb += std::log(dgt.genome.ref_pprob);

    /// TODO STREL-125 find a way to restore strand bias feature
    // allele.strandBias=dgt.strand_bias;

    updateSiteSampleInfo(opt, sampleIndex, cpi, locus);
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

    // prep step 1) clean pileups in all samples:
    for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
    {
        _pileupCleaner.CleanPileupErrorProb(sample(sampleIndex).cpi);
    }

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
        const int ploidy(get_ploidy(pos, sampleIndex));
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
    getSiteAltAlleles(refBaseIndex, allDgt, altAlleles);

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
    for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
    {
        updateSnvLocusWithSampleInfo(
            _opt, sample(sampleIndex), callerPloidy[sampleIndex], groupLocusPloidy[sampleIndex],
            allDgt[sampleIndex], sampleIndex, *locusPtr, homRefLogProb);
    }

    // add sample-independent info:
    locusPtr->anyVariantAlleleQuality = ln_error_prob_to_qphred(homRefLogProb);

    if (isForcedOutput or locusPtr->isVariantLocus())
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
    const unsigned baseId, ///< TODO STREL-125 TMP!!!!!!
    GermlineContinuousSiteLocusInfo& locus)
{
    const auto& siteSampleInfo(locus.getSiteSample(sampleIndex));

    GermlineContinuousSiteSampleInfo siteContinuousSampleInfo;

    siteContinuousSampleInfo.continuousTotalDepth = siteSampleInfo.spanningDeletionReadCount;
    for (unsigned baseId2(0); baseId2 < N_BASE; ++baseId2)
    {
        /// TODO STREL-125 sample generalization
        siteContinuousSampleInfo.continuousTotalDepth += siteSampleInfo.alleleObservationCounts(baseId2);
    }

    siteContinuousSampleInfo.continuousAlleleDepth = siteSampleInfo.alleleObservationCounts(baseId);

    locus.setContinuousSiteSampleInfo(sampleIndex, siteContinuousSampleInfo);
}



/// strand bias values summed over all samples
struct StrandBiasCounts
{
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
    const PileupCleaner& pileupCleaner,
    const unsigned sampleIndex,
    const unsigned baseIndex, ///< TODO STREL-125 TMP!!!!!!
    std::vector<StrandBiasCounts>& strandBiasCounts,
    GermlineContinuousSiteLocusInfo& locus)
{
    const bool isRefAllele(locus.refBaseIndex == baseIndex);

    auto& siteAlleles(locus.getSiteAlleles());
    const auto alleleCount(siteAlleles.size());

    const CleanedPileup& cpi(sif.cpi);

    updateSiteSampleInfo(opt, sampleIndex, cpi, locus);
    updateContinuousSiteSampleInfo(sampleIndex, baseIndex, locus);


    //const snp_pos_info& pi(cpi.rawPileup());

    pileupCleaner.CleanPileupErrorProb(sif.cpi);

    const snp_pos_info& good_pi(cpi.cleanedPileup());
//    const extended_pos_info& good_epi(cpi.getExtendedPosInfo());

    {
        auto& sampleInfo(locus.getSample(sampleIndex));
        const auto& continuousSiteSampleInfo(locus.getContinuousSiteSample(sampleIndex));

        const double alleleFrequency(continuousSiteSampleInfo.getContinuousAlleleFrequency());

        {
            // use diploid gt codes as a convenient way to summarize the continuous variant calls:
            static const VcfGenotype homrefGtIndex(0, 0);
            static const VcfGenotype hetGtIndex(0, 1);
            static const VcfGenotype homGtIndex(1, 1);

            auto getGtIndex = [&]() -> VcfGenotype {
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
                (unsigned) opt.min_qscore, 40);

        // update strand bias intermediates:
        for (unsigned alleleIndex(0); alleleIndex < alleleCount; ++alleleIndex)
        {
            const auto& allele(siteAlleles[alleleIndex]);
            auto& sbcounts(strandBiasCounts[alleleIndex]);

            for (const base_call& bc : good_pi.calls)
            {
                if (bc.is_fwd_strand)
                {
                    if (bc.base_id == allele.baseIndex)
                        sbcounts.fwdAlt++;
                    else
                        sbcounts.fwdOther++;
                }
                else if (bc.base_id == allele.baseIndex)
                    sbcounts.revAlt++;
                else
                    sbcounts.revOther++;
            }
        }
    }
}



void
updateContinuousSnvLocusInfo(
    const starling_options& opt,
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

    // update strand bias for each allele:
    for (unsigned alleleIndex(0); alleleIndex < alleleCount; ++alleleIndex)
    {
        auto& allele(siteAlleles[alleleIndex]);
        const auto& sbcounts(strandBiasCounts[alleleIndex]);
        allele.strandBias = starling_continuous_variant_caller::strand_bias(
            sbcounts.fwdAlt, sbcounts.revAlt, sbcounts.fwdOther, sbcounts.revOther, opt.noise_floor);
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

    // end sample generalization
    /// TODO STREL-125 generalize to multisample
    assert(sampleCount == 1);
    const unsigned tmpSampleIndex(0);

    sample_info& sif(sample(tmpSampleIndex));

    const CleanedPileup& cpi(sif.cpi);
    const snp_pos_info& pi(cpi.rawPileup());

    _pileupCleaner.CleanPileupErrorProb(sif.cpi);

    const uint8_t refBaseId = base_to_id(pi.get_ref_base());

    // report one locus (ie. vcf record) per alt allele in continuous mode
    bool isAnySiteOutputAtPosition(false);

    auto addBase = [&](const uint8_t baseId, const bool isForcedOutputUsed)
    {
        const bool isRefAllele(baseId == refBaseId);
        std::unique_ptr<GermlineContinuousSiteLocusInfo> locusPtr(new GermlineContinuousSiteLocusInfo(
            sampleCount, pos, base_to_id(pi.get_ref_base()), isForcedOutput));

        // setup alt allele first:
        if (not isRefAllele)
        {
            locusPtr->addAltSiteAllele(static_cast<BASE_ID::index_t>(baseId));
        }

        // set some sample-dependent info
        const auto& siteAlleles(locusPtr->getSiteAlleles());
        const auto alleleCount(siteAlleles.size());
        std::vector<StrandBiasCounts> strandBiasCounts(alleleCount);

        for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
        {
            updateContinuousSnvLocusWithSampleInfo(
                _opt, sample(sampleIndex), _pileupCleaner, sampleIndex, baseId, strandBiasCounts, *locusPtr);
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

    for (unsigned baseId(0); baseId < N_BASE; ++baseId)
    {
        if (baseId == refBaseId) continue;
        addBase(baseId,isForcedOutput);
    }

    // ensure that at least one base is added for site
    if (not isAnySiteOutputAtPosition)
    {
        addBase(refBaseId,true);
    }
}



void
starling_pos_processor::
process_pos_indel(const pos_t pos)
{
    if (_opt.is_bsnp_diploid())
    {
        process_pos_indel_digt(pos);
    }
    else
    {
        process_pos_indel_continuous(pos);
    }
}



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



/// setup indelSampleInfo assuming that corresponding sampleInfo has already been initialized
void
updateIndelSampleInfo(
    const starling_base_options& opt,
    const starling_deriv_options& dopt,
    const pos_basecall_buffer& basecallBuffer,
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    const unsigned sampleIndex,
    GermlineIndelLocusInfo& locus)
{
    GermlineIndelSampleInfo indelSampleInfo;

    const LocusSampleInfo& sampleInfo(locus.getSample(sampleIndex));
    const unsigned callerPloidy(sampleInfo.getPloidy().getPloidy());

    // set sitePloidy:
    const auto& range(locus.range());

    auto& sitePloidy(indelSampleInfo.sitePloidy);
    sitePloidy.resize(range.size());

    auto updateSitePloidyForAlleleIndex = [&](const uint8_t alleleIndex) {
        if (alleleIndex == 0) return;

        const IndelKey& indelKey2(locus.getIndelAlleles()[alleleIndex - 1].indelKey);

        pos_t leadingOffset(indelKey2.pos - range.begin_pos());
        pos_t trailingOffset(indelKey2.right_pos() - range.begin_pos());
        for (pos_t locusOffset(leadingOffset); locusOffset < trailingOffset; ++locusOffset)
        {
            sitePloidy[locusOffset] -= 1;
        }
    };

    std::fill(sitePloidy.begin(), sitePloidy.end(), callerPloidy);
    const auto& maxGt(sampleInfo.max_gt());
    if (callerPloidy == 2)
    {
        updateSitePloidyForAlleleIndex(maxGt.getAllele0Index());
        updateSitePloidyForAlleleIndex(maxGt.getAllele0Index());
    }
    else if (callerPloidy == 1)
    {
        updateSitePloidyForAlleleIndex(maxGt.getAllele0Index());
    }
    else
    {
        assert(false);
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
static
void
updateIndelLocusWithSampleInfo(
    const starling_base_options& opt,
    const starling_deriv_options& dopt,
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    const unsigned topAlleleIndexInSample,
    const OrthogonalVariantAlleleCandidateGroup& contrastGroup,
    const starling_pos_processor::sample_info& sif,
    const unsigned callerPloidy,
    const unsigned groupLocusPloidy,
    const unsigned sampleIndex,
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
        opt, dopt, sif.sample_opt, callerPloidy, sampleIndex, alleleGroup, contrastGroup,
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
            // don't enforce maxQ at this point, b/c we're going to possibly select down from this list:
            samplePL[gtIndex] = ln_error_prob_to_qphred(genotypeLogLhood[gtIndex]-genotypeLogLhood[maxGtIndex]);
        }
    }

    //------------------------------------------------
    // compute posteriors/qualities:

    // get patternRepeatCount, if more than one alt allele then base this off of the most likely allele per sample
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
        for (unsigned allele0Index(0); allele0Index<=fullAlleleCount; ++allele0Index)
        {
            const unsigned genotypeIndex(VcfGenotypeUtil::getGenotypeIndex(allele0Index));
            const unsigned priorIndex(getPriorIndex(topAlleleIndexInSample, allele0Index));
            genotypePosterior[genotypeIndex] = genotypeLogLhood[genotypeIndex] + genotypeLogPrior[priorIndex];
        }
    }
    else if(callerPloidy == 2)
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

    setGentypeFromGenotypeIndex(callerPloidy, maxGenotypeIndex, sampleInfo.maxGenotypeIndexPolymorphic);

    // get genome posterior
    if (callerPloidy == 1)
    {
        static const bool isHaploid(true);
        const double* genotypeLogPrior(genotypePriors.getNAllele(isHaploid));
        for (unsigned allele0Index(0); allele0Index<=fullAlleleCount; ++allele0Index)
        {
            const unsigned genotypeIndex(VcfGenotypeUtil::getGenotypeIndex(allele0Index));
            const unsigned priorIndex(getPriorIndex(topAlleleIndexInSample, allele0Index));
            genotypePosterior[genotypeIndex] = genotypeLogLhood[genotypeIndex] + genotypeLogPrior[priorIndex];
        }
    }
    else if(callerPloidy == 2)
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

    setGentypeFromGenotypeIndex(callerPloidy, maxGenotypeIndex, sampleInfo.maxGenotypeIndex);

    // set GQX
    sampleInfo.setGqx();

    // update homref prob for QUAL
    homRefLogProb += std::log(genotypePosterior[AG_GENOTYPE::HOMREF]);

    // set indelSampleInfo
    updateIndelSampleInfo(opt, dopt, sif.bc_buff, alleleGroup, sampleIndex, locus);
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

    // track all forced output alleles in a separate group (even if they go into topVariant group)
    // to ensure that these are output even if not included in the most likely genotype for any sample
    //
    OrthogonalVariantAlleleCandidateGroup forcedOutputAlleleGroup;
    {
        const unsigned orthogonalVariantAlleleCount(orthogonalVariantAlleles.size());
        for (unsigned alleleIndex(0); alleleIndex < orthogonalVariantAlleleCount; alleleIndex++)
        {
            const IndelData& indelData(orthogonalVariantAlleles.data(alleleIndex));
            if (indelData.isForcedOutput)
            {
                forcedOutputAlleleGroup.addVariantAllele(orthogonalVariantAlleles.iter(alleleIndex));
            }
        }
    }

    // track top alt allele within each sample -- this is used as a temporary crutch to transfer the previous prior
    // calcualation from single to multi-sample, and should be removed when a matured prior scheme is put in place
    std::vector<unsigned> topVariantAlleleIndexPerSample(sampleCount);

    // rank input alleles to pick the top N, N=ploidy, per sample, and aggregate/rank these
    // over all samples
    OrthogonalVariantAlleleCandidateGroup topVariantAlleleGroup;
    selectTopOrthogonalAllelesInAllSamples(
        sampleCount, callerPloidy, orthogonalVariantAlleles, topVariantAlleleGroup, topVariantAlleleIndexPerSample);

    // At this point topVariantAlleleGroup represents the best alleles which
    // start at the current position (over all samples). Now we add conflicting
    // alleles at other positions and re-rank, re-select the top alleles again:
    if (not topVariantAlleleGroup.empty())
    {
        addAllelesAtOtherPositions(sampleCount, callerPloidy, pos, get_largest_total_indel_ref_span_per_read(),
                                   getIndelBuffer(), topVariantAlleleGroup, topVariantAlleleIndexPerSample);
    }

    // genotype and report topVariantAlleleGroup
    //

    // overlapping allele groups are reported only once, when grouped together from the left-most position
    bool isReportableLocus(true);
    bool isReportedLocus(false);
    {
        const unsigned alleleGroupSize(topVariantAlleleGroup.size());
        for (unsigned genotypeAlleleIndex(0); genotypeAlleleIndex < alleleGroupSize; ++genotypeAlleleIndex)
        {
            const IndelKey& indelKey(topVariantAlleleGroup.key(genotypeAlleleIndex));
            if (indelKey.pos < pos)
            {
                isReportableLocus=false;
                break;
            }
        }
    }

    if (isReportableLocus)
    {
        // parameter inputs if/when we wrap this as a function:
        static const bool isForcedOutput(false);
        static OrthogonalVariantAlleleCandidateGroup emptyGroup;

        // setup new indel locus:
        std::unique_ptr<GermlineIndelLocusInfo> locusPtr(new GermlineDiploidIndelLocusInfo(_dopt.gvcf, sampleCount));

        // cycle through variant alleles and add them to locus (the locus interface requires that this is done first):
        addIndelAllelesToLocus(topVariantAlleleGroup, isForcedOutput, *locusPtr);

        // add sample-dependent info:
        double homRefLogProb(0);
        for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
        {
            updateIndelLocusWithSampleInfo(
                _opt, _dopt, topVariantAlleleGroup, topVariantAlleleIndexPerSample[sampleIndex], emptyGroup, sample(sampleIndex),
                callerPloidy[sampleIndex], groupLocusPloidy[sampleIndex], sampleIndex, *locusPtr, homRefLogProb);
        }

        // add sample-independent info:
        locusPtr->anyVariantAlleleQuality = ln_error_prob_to_qphred(homRefLogProb);

        if (isForcedOutput or locusPtr->isVariantLocus())
        {
            // finished! send this locus down the pipe:
            _gvcfer->add_indel(std::move(locusPtr));

            isReportedLocus = true;
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
        // trim the forced output allele set to take out any alleles already called as variants:
        if (not forcedOutputAlleleGroup.empty())
        {
            const unsigned forcedCount(forcedOutputAlleleGroup.size());
            for (unsigned forcedIndex(0); forcedIndex < forcedCount; ++forcedIndex)
            {
                const unsigned reverseForcedIndex(forcedCount-(forcedIndex+1));
                if (_forcedAllelesAlreadyOutput.count(forcedOutputAlleleGroup.key(reverseForcedIndex)) > 0)
                {
                    forcedOutputAlleleGroup.alleles.erase(forcedOutputAlleleGroup.alleles.begin() + reverseForcedIndex);
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
            for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
            {
                updateIndelLocusWithSampleInfo(
                    _opt, _dopt, fakeForcedOutputAlleleGroup, fakeTopVariantAlleleIndexPerSample, topVariantAlleleGroup, sample(sampleIndex),
                    callerPloidy[sampleIndex], groupLocusPloidy[sampleIndex], sampleIndex, *locusPtr, homRefLogProb);
            }

            // add sample-independent info:
            locusPtr->anyVariantAlleleQuality = ln_error_prob_to_qphred(homRefLogProb);

            // finished! send this locus down the pipe:
            _gvcfer->add_indel(std::move(locusPtr));
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
            (unsigned) opt.min_qscore, 40);

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



void
starling_pos_processor::
write_counts(const pos_range& output_report_range) const
{
    std::ostream* report_osptr(get_report_osptr());
    if (NULL==report_osptr) return;
    std::ostream& report_os(*report_osptr);

    const sample_info& sif(sample());

    report_os << std::setprecision(8);
    report_stream_stat(sif.ss,"ALLSITES_COVERAGE",output_report_range,report_os);
    report_stream_stat(sif.used_ss,"ALLSITES_COVERAGE_USED",output_report_range,report_os);

    if (_opt.is_ref_set())
    {
        report_stream_stat(sif.ssn,"NO_REF_N_COVERAGE",output_report_range,report_os);
        report_stream_stat(sif.used_ssn,"NO_REF_N_COVERAGE_USED",output_report_range,report_os);
    }
}

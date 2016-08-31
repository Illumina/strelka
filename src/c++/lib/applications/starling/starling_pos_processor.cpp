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



static
void
report_counts(
    const snp_pos_info& pi,
    const unsigned n_unused_calls,
    const pos_t output_pos,
    std::ostream& os)
{
    unsigned base_count[N_BASE];

    for (unsigned i(0); i<N_BASE; ++i) base_count[i] = 0;

    for (const auto& bc : pi.calls)
    {
        assert(bc.base_id!=BASE_ID::ANY);
        base_count[bc.base_id]++;
    }

    os << output_pos << '\t';
    for (unsigned i(0); i<N_BASE; ++i)
    {
        os << base_count[i] << '\t';
    }
    os << n_unused_calls << '\n';
}



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



void
starling_pos_processor::
process_pos_snp_continuous(const pos_t pos)
{
    const unsigned sampleCount(getSampleCount());
    const bool isForcedOutput(is_forced_output_pos(pos));

    // end sample generalization
    /// TODO STREL-125 generalize to multisample
    assert(sampleCount == 1);
    const unsigned sampleIndex(0);

    sample_info& sif(sample(sampleIndex));

    const CleanedPileup& cpi(sif.cpi);
    const snp_pos_info& pi(cpi.rawPileup());

    _pileupCleaner.CleanPileupErrorProb(sif.cpi);

    const snp_pos_info& good_pi(cpi.cleanedPileup());

    GermlineContinuousSiteLocusInfo templateLocus(sampleCount, pos, pi.get_ref_base(), good_pi,
                                            _opt.used_allele_count_min_qscore, _opt.min_het_vf, isForcedOutput);

    templateLocus.n_used_calls = cpi.n_used_calls();
    templateLocus.n_unused_calls = cpi.n_unused_calls();
    // hpol filter
    templateLocus.hpol = get_snp_hpol_size(pos, _ref);

    if (_opt.is_counts)
    {
        report_counts(good_pi, templateLocus.n_unused_calls, templateLocus.pos + 1, *_streams.counts_osptr());
    }

    // report one locus (ie. vcf record) per alt allele in continuous mode
    bool isSiteAddedForPosition(false);

    auto addBase = [&](const uint8_t baseId, const bool isForcedOutputUsed)
    {
        std::unique_ptr<GermlineContinuousSiteLocusInfo> locusPtr(new GermlineContinuousSiteLocusInfo(templateLocus));
        starling_continuous_variant_caller::position_snp_call_continuous(_opt, good_pi, baseId, isForcedOutputUsed,
                                                                         (GermlineContinuousSiteLocusInfo&) *locusPtr);
        if (not locusPtr->altAlleles.empty())
        {
            isSiteAddedForPosition = true;
//            _gvcfer->add_site(std::move(si));
        }
    };

    for (unsigned baseId(0); baseId < N_BASE; ++baseId)
    {
        addBase(baseId,isForcedOutput);
    }

    /// ensure that at least one base is added for site
    if (not isSiteAddedForPosition)
    {
        const uint8_t refBaseId = base_to_id(templateLocus.ref);
        addBase(refBaseId,true);
    }
}



void
starling_pos_processor::
process_pos_snp_digt(const pos_t pos)
{
    const unsigned sampleCount(getSampleCount());
    const bool isForcedOutput(is_forced_output_pos(pos));

    /// end multi-sample generalization:
    assert(sampleCount == 1);
    const unsigned sampleIndex(0);
    sample_info& sif(sample(sampleIndex));

    const CleanedPileup& cpi(sif.cpi);
    const snp_pos_info& pi(cpi.rawPileup());

    _pileupCleaner.CleanPileupErrorProb(sif.cpi);

    const snp_pos_info& good_pi(cpi.cleanedPileup());
    const extended_pos_info& good_epi(cpi.getExtendedPosInfo());

    std::unique_ptr<GermlineDiploidSiteLocusInfo> locusPtr(new GermlineDiploidSiteLocusInfo(_dopt.gvcf, sampleCount, pos,pi.get_ref_base(),good_pi,_opt.used_allele_count_min_qscore, isForcedOutput));
    locusPtr->n_used_calls=cpi.n_used_calls();
    locusPtr->n_unused_calls=cpi.n_unused_calls();

    // check whether we're in a haploid region:
    locusPtr->dgt.ploidy=(get_ploidy(pos, sampleIndex));

    const pos_t output_pos(pos+1);

    if (_opt.is_counts)
    {
        report_counts(good_pi,locusPtr->n_unused_calls,output_pos,*_streams.counts_osptr());
    }

    if (_opt.is_bsnp_diploid())
    {
        _dopt.pdcaller().position_snp_call_pprob_digt(
            _opt,good_epi,locusPtr->dgt, _opt.is_all_sites());
    }

    //    const bool is_nf_snp(is_snp && (! is_filter_snp));
    if (isForcedOutput or locusPtr->dgt.is_snp())
    {
        if (_opt.is_compute_hapscore)
        {
            locusPtr->hapscore=get_hapscore(pi.hap_set);
        }

        // calculate empirical scoring metrics
        if (_opt.is_compute_germline_scoring_metrics())
        {
            locusPtr->mapqRMS = pi.mapqTracker.getRMS();
            locusPtr->mapqZeroCount = pi.mapqTracker.zeroCount;
            locusPtr->mapqCount = pi.mapqTracker.count;
            locusPtr->ReadPosRankSum = pi.get_read_pos_ranksum();
            locusPtr->MQRankSum = pi.get_mq_ranksum();
            locusPtr->BaseQRankSum = pi.get_baseq_ranksum();
            locusPtr->rawPos = pi.get_raw_pos();
            locusPtr->avgBaseQ = pi.get_raw_baseQ();
        }

        // hpol filter
        locusPtr->hpol = get_snp_hpol_size(pos,_ref);
    }

    //Add site to gvcf
//    _gvcfer->add_site(std::move(si));
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
    if (callerPloidy == 2)
    {
        uint8_t allele0Index, allele1Index;
        VcfGenotypeUtil::getAlleleIndices(sampleInfo.maxGenotypeIndexPolymorphic, allele0Index,
                                          allele1Index);
        updateSitePloidyForAlleleIndex(allele0Index);
        updateSitePloidyForAlleleIndex(allele1Index);
    }
    else if (callerPloidy == 1)
    {
        uint8_t allele0Index;
        VcfGenotypeUtil::getAlleleIndices(sampleInfo.maxGenotypeIndexPolymorphic, allele0Index);
        updateSitePloidyForAlleleIndex(allele0Index);
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

    normalize_ln_distro(std::begin(genotypePosterior), std::end(genotypePosterior), sampleInfo.maxGenotypeIndexPolymorphic);

    sampleInfo.genotypeQualityPolymorphic = error_prob_to_qphred(prob_comp(std::begin(genotypePosterior), std::end(genotypePosterior), sampleInfo.maxGenotypeIndexPolymorphic));


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

    normalize_ln_distro(std::begin(genotypePosterior), std::end(genotypePosterior), sampleInfo.maxGenotypeIndex);

    sampleInfo.genotypeQuality = error_prob_to_qphred(prob_comp(std::begin(genotypePosterior), std::end(genotypePosterior), sampleInfo.maxGenotypeIndex));

    // set GQX
    // maxGenotypeIndex != maxGenotypeIndexPolymorphic indicates we're in a boundary zone between variant and hom-ref call
    if (sampleInfo.maxGenotypeIndex != sampleInfo.maxGenotypeIndexPolymorphic)
    {
        sampleInfo.gqx = 0;
    }
    else
    {
        sampleInfo.gqx = std::min(sampleInfo.genotypeQuality, sampleInfo.genotypeQualityPolymorphic);
    }

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

        if (isForcedOutput or (locusPtr->anyVariantAlleleQuality != 0))
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
    static const unsigned hetGtIndex(VcfGenotypeUtil::getGenotypeIndex(0,1));
    static const unsigned homGtIndex(VcfGenotypeUtil::getGenotypeIndex(1,1));

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

        // determine if this locus is printable:
        bool isReportableAllele(isForcedOutput);
        for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
        {
            const auto& indelSampleInfo(locusPtr->getIndelSample(sampleIndex));
            if (not isReportableAllele)
            {
                const double alleleFrequency(indelSampleInfo.alleleFrequency());
                if ((indelSampleInfo.legacyReportInfo.n_confident_indel_reads > 0) and
                    (alleleFrequency > _opt.min_het_vf))
                {
                    isReportableAllele = true;
                }
            }
        }
        if (not isReportableAllele) continue;

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

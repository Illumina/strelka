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
#include "blt_common/position_nonref_test.hh"
#include "blt_common/position_nonref_2allele_test.hh"
#include "blt_common/ref_context.hh"
#include "blt_util/log.hh"
#include "starling_common/AlleleGroupGenotype.hh"
#include "starling_common/AlleleReportInfoUtil.hh"
#include "starling_common/indel_util.hh"
#include "starling_common/OrthogonalVariantAlleleCandidateGroup.hh"
#include "starling_common/OrthogonalVariantAlleleCandidateGroupUtil.hh"

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



static
void
write_snp_prefix_info_file(const std::string& seq_name,
                           const pos_t output_pos,
                           const char ref,
                           const unsigned n_used_calls,
                           const unsigned n_unused_calls,
                           std::ostream& os)
{
    os << seq_name << "\t"
       << output_pos << "\t"
       << n_used_calls << "\t"
       << n_unused_calls << "\t"
       << ref;
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
process_pos_snp_single_sample(
    const pos_t pos,
    const unsigned sample_no)
{
    try
    {
        if (_opt.is_bsnp_diploid())
        {
            process_pos_snp_single_sample_impl(pos,sample_no);
        }
        else
        {
            process_pos_snp_single_sample_continuous(pos,sample_no);
        }

    }
    catch (...)
    {
        log_os << "Exception caught in starling_pos_processor_base.process_pos_snp_single_sample_impl() while processing chromosome position: " << (pos+1) << "\n"
               << "snp_pos_info:\n";
        log_os << sample(sample_no).bc_buff.get_pos(pos) << "\n";
        throw;
    }
}

void
starling_pos_processor::
process_pos_snp_single_sample_continuous(
    const pos_t pos,
    const unsigned sample_no)
{
    assert(_opt.gvcf.is_gvcf_output());

    if (sample_no!=0) return;

    /// TODO STREL-125 generalize to multisample
    const unsigned sampleCount(1);

    sample_info& sif(sample(sample_no));

    const CleanedPileup& cpi(sif.cpi);
    const snp_pos_info& pi(cpi.rawPileup());

    _pileupCleaner.CleanPileupErrorProb(sif.cpi);

    const snp_pos_info& good_pi(cpi.cleanedPileup());
    const bool is_forced(is_forced_output_pos(pos));

    if (pi.calls.empty() && !is_forced) return;

    GermlineContinuousSiteLocusInfo locusInfo(sampleCount, pos, pi.get_ref_base(), good_pi,
                                            _opt.used_allele_count_min_qscore, _opt.min_het_vf, is_forced);

    locusInfo.n_used_calls = cpi.n_used_calls();
    locusInfo.n_unused_calls = cpi.n_unused_calls();
    // hpol filter
    locusInfo.hpol = get_snp_hpol_size(pos, _ref);

    if (_opt.is_counts)
    {
        report_counts(good_pi, locusInfo.n_unused_calls, locusInfo.pos + 1, *_streams.counts_osptr());
    }

    // report one locus (ie. vcf record) per alt allele in continuous mode
    bool isSiteAddedForPosition(false);

    auto addBase = [&](const uint8_t baseId, const bool isForcedOutput)
    {
        std::unique_ptr<GermlineContinuousSiteLocusInfo> si(new GermlineContinuousSiteLocusInfo(locusInfo));
        starling_continuous_variant_caller::position_snp_call_continuous(_opt, good_pi, baseId, isForcedOutput,
                                                                         (GermlineContinuousSiteLocusInfo&) *si);
        if (not si->altAlleles.empty())
        {
            isSiteAddedForPosition = true;
            _gvcfer->add_site(std::move(si));
        }
    };

    for (unsigned baseId(0); baseId < N_BASE; ++baseId)
    {
        addBase(baseId,is_forced);
    }

    /// ensure that at least one base is added for site
    if (not isSiteAddedForPosition)
    {
        const uint8_t refBaseId = base_to_id(locusInfo.ref);
        addBase(refBaseId,true);
    }
}



void
starling_pos_processor::
process_pos_snp_single_sample_impl(
    const pos_t pos,
    const unsigned sampleIndex)
{
    // TODO:
    //
    // note this might not matter wrt larger changes taking place, but here goes:
    //
    // change filters to support vcf concept of 1..N filters which are added to the genotype information
    //
    // generalize site tests with an object
    //
    // genotype_test {
    //    ctor(); // setup any cached values
    //
    //    test(site_info);
    //
    //    write()?? (do we need to even bother with this?)
    // }
    //

    sample_info& sif(sample(sampleIndex));

    const CleanedPileup& cpi(sif.cpi);
    const snp_pos_info& pi(cpi.rawPileup());

    // note multi-sample status -- can still be called only for one sample
    // and only for sample 0. working on generalization:
    //
    if (sampleIndex!=0) return;

    const unsigned sampleCount(1);

    const bool is_forced(is_forced_output_pos(pos));

    // the second term in is_skippable below forces sites to go through the pipeline
    // if phaser has put a hold on buffer cleanup. This ensures that the phaser will be turned back off
    //
    // TODO: there must be a way to force correct usage into the phaser's API instead of requiring this brittle hack
    const bool is_skippable(! (is_forced || is_save_pileup_buffer()));

    if (pi.calls.empty() && is_skippable) return;

    _pileupCleaner.CleanPileupErrorProb(sif.cpi);

    const snp_pos_info& good_pi(cpi.cleanedPileup());
    const extended_pos_info& good_epi(cpi.getExtendedPosInfo());


    std::unique_ptr<GermlineDiploidSiteLocusInfo> si(new GermlineDiploidSiteLocusInfo(_dopt.gvcf, sampleCount, pos,pi.get_ref_base(),good_pi,_opt.used_allele_count_min_qscore, is_forced));
    si->n_used_calls=cpi.n_used_calls();
    si->n_unused_calls=cpi.n_unused_calls();


    // delay writing any snpcalls so that anomaly tests can (optionally) be applied as filters:
    //
    nonref_test_call nrc;
    //lrt_snp_call lsc;
    //std::unique_ptr<nploid_genotype> ngt_ptr;

    // check whether we're in a haploid region:
    si->dgt.ploidy=(get_ploidy(pos, sampleIndex));

    const pos_t output_pos(pos+1);

    if (_opt.is_counts)
    {
        report_counts(good_pi,si->n_unused_calls,output_pos,*_streams.counts_osptr());
    }

    if (_opt.is_nonref_test() || _opt.is_nonref_sites())
    {
        position_nonref_2allele_test(good_pi,
                                     _opt,
                                     _opt.is_nonref_sites(),
                                     nrc);
#if 0
        static const bool is_mle_freq(false);

        position_nonref_test(good_pi,
                             _opt.nonref_variant_rate,
                             _opt.min_nonref_freq,
                             is_mle_freq,
                             nrc);
#endif

    }

#if 0
    if (_opt.is_lsnp)
    {
        position_snp_call_lrt(_opt.lsnp_alpha,good_pi,lsc);
    }
#endif
    if (_opt.is_bsnp_diploid())
    {
        _dopt.pdcaller().position_snp_call_pprob_digt(
            _opt,good_epi,si->dgt, _opt.is_all_sites());
    }
#if 0
    if (_opt.is_bsnp_nploid)
    {
        ngt_ptr.reset(new nploid_genotype(*_ninfo));
        position_snp_call_pprob_nploid(_opt.bsnp_nploid_snp_prob,good_pi,*_ninfo,*ngt_ptr);
    }
#endif

    const bool is_snp(nrc.is_snp || si->dgt.is_snp);

    //    const bool is_nf_snp(is_snp && (! is_filter_snp));
    if (is_snp || is_forced)
    {
        if (_opt.is_compute_hapscore)
        {
            si->hapscore=get_hapscore(pi.hap_set);
        }

        // calculate empirical scoring metrics
        if (_opt.is_compute_germline_scoring_metrics())
        {
            si->mapqRMS = pi.mapqTracker.getRMS();
            si->mapqZeroCount = pi.mapqTracker.zeroCount;
            si->mapqCount = pi.mapqTracker.count;
            si->ReadPosRankSum = pi.get_read_pos_ranksum();
            si->MQRankSum = pi.get_mq_ranksum();
            si->BaseQRankSum = pi.get_baseq_ranksum();
            si->rawPos = pi.get_raw_pos();
            si->avgBaseQ = pi.get_raw_baseQ();
        }

        // hpol filter
        si->hpol = get_snp_hpol_size(pos,_ref);
    }

    if (_opt.is_nonref_sites())
    {
        std::ostream& bos(*_streams.nonref_sites_osptr());
        write_snp_prefix_info_file(_chrom_name,output_pos,pi.get_ref_base(),si->n_used_calls,si->n_unused_calls,bos);
        bos << "\t";
        write_nonref_2allele_test(_opt,good_pi,nrc,bos);
        bos << "\n";
    }

    // report events:
    //
    bool is_reported_event(false);

    std::ostream& report_os(std::cerr);

    if (is_snp)
    {
        if (nrc.is_snp)
        {
            std::ostream& bos(*_streams.nonref_test_osptr());
            write_snp_prefix_info_file(_chrom_name,output_pos,pi.get_ref_base(),si->n_used_calls,si->n_unused_calls,bos);
            bos << "\t";
            write_nonref_2allele_test(_opt,good_pi,nrc,bos);
#if 0
            write_nonref_test(_opt,good_pi,nrc,bos);
#endif
            bos << "\n";
        }
#if 0
        if (lsc.is_snp)
        {
            write_snp_prefix_info("LSNP",output_pos,pi.ref_base,_site_info.n_used_calls,_site_info.n_unused_calls,report_os);
            report_os << " " << lsc << "\n";
        }
        if (mgt.is_snp)
        {
            write_snp_prefix_info("BSNP1",output_pos,pi.ref_base,_site_info.n_used_calls,_site_info.n_unused_calls,report_os);
            report_os << " " << mgt << "\n";
        }
        if (ngt_ptr.get() && ngt_ptr->is_snp)
        {
            write_snp_prefix_info("BSNPN",output_pos,pi.ref_base,_site_info.n_used_calls,_site_info.n_unused_calls,report_os);
            report_os << " ";
            nploid_write(*_ninfo,*ngt_ptr,report_os);
            report_os << "\n";
        }
#endif

        is_reported_event = true;
    }

    if (_opt.is_print_all_site_evidence || (_opt.is_print_evidence && is_reported_event))
    {
        report_os << "EVIDENCE pos: " << output_pos << "\n"
                  << "is_snp: " << is_snp << "\n"
                  << pi << "\n";
    }

    //Add site to gvcf
    if (_opt.gvcf.is_gvcf_output())
    {
        _gvcfer->add_site(std::move(si));
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



/// translate AG_GENOTYPE to older genotype indices
///
/// TODO retire the old STAR_DIINDEL scheme and remove these conversion routines
namespace AG_GENOTYPE
{
static
STAR_DIINDEL::index_t
mapAlleleToDindel(
    const unsigned genotypeId,
    const unsigned alleleId)
{
    if (alleleId == 0)
    {
        switch (static_cast<index_t>(genotypeId))
        {
        case HOM0:
            return STAR_DIINDEL::HOM;
        case HET0:
        case HET01:
            return STAR_DIINDEL::HET;
        default:
            return STAR_DIINDEL::NOINDEL;
        }
    }
    else if (alleleId == 1)
    {
        switch (static_cast<index_t>(genotypeId))
        {
        case HOM1:
            return STAR_DIINDEL::HOM;
        case HET1:
        case HET01:
            return STAR_DIINDEL::HET;
        default:
            return STAR_DIINDEL::NOINDEL;
        }
    }
    else
    {
        assert(false and "Unknown alleleId");
        return STAR_DIINDEL::SIZE;
    }
}
}



/// TODO remove this function once we eliminate starling_diploid_indel as an intermediary format
static
starling_diploid_indel
locusGenotypeToDindel(
    const AlleleGroupGenotype& locusGenotype,
    const unsigned alleleId)
{
    starling_diploid_indel dindel;

    if (locusGenotype.maxGenotypeIndex == AG_GENOTYPE::HET01)
    {
        dindel.is_diplotype_model_hetalt = true;
    }

    dindel.pprob[STAR_DIINDEL::NOINDEL] = locusGenotype.posteriorProb[AG_GENOTYPE::HOMREF];
    dindel.pprob[STAR_DIINDEL::HET] = locusGenotype.posteriorProb[AG_GENOTYPE::getAlleleHetId(alleleId)];
    dindel.pprob[STAR_DIINDEL::HOM] = locusGenotype.posteriorProb[AG_GENOTYPE::getAlleleHomId(alleleId)];

    if (locusGenotype.maxGenotypeIndex == AG_GENOTYPE::HET01)
    {
        dindel.pprob[STAR_DIINDEL::HET] = locusGenotype.posteriorProb[AG_GENOTYPE::HET01];
    }

    dindel.max_gt = AG_GENOTYPE::mapAlleleToDindel(locusGenotype.maxGenotypeIndex, alleleId);
    dindel.max_gt_qphred = (int) locusGenotype.genotypeQuality;
    dindel.max_gt_poly = AG_GENOTYPE::mapAlleleToDindel(locusGenotype.maxGenotypeIndexPolymorphic, alleleId);
    dindel.max_gt_poly_qphred = (int) locusGenotype.genotypeQualityPolymorphic;

    dindel.phredLoghood[STAR_DIINDEL::NOINDEL] = (unsigned) locusGenotype.phredLoghood[AG_GENOTYPE::HOMREF];
    dindel.phredLoghood[STAR_DIINDEL::HET] = (unsigned) locusGenotype.phredLoghood[AG_GENOTYPE::getAlleleHetId(alleleId)];
    dindel.phredLoghood[STAR_DIINDEL::HOM] = (unsigned) locusGenotype.phredLoghood[AG_GENOTYPE::getAlleleHomId(alleleId)];

    if (locusGenotype.maxGenotypeIndex == AG_GENOTYPE::HET01)
    {
        dindel.phredLoghood[STAR_DIINDEL::HET] = (unsigned) locusGenotype.phredLoghood[AG_GENOTYPE::HET01];
    }

    {
        unsigned minIndex(0);
        for (unsigned gt(1); gt < STAR_DIINDEL::SIZE; ++gt)
        {
            if (dindel.phredLoghood[gt] < dindel.phredLoghood[minIndex]) minIndex = gt;
        }
        const int minPL(dindel.phredLoghood[minIndex]);
        for (unsigned gt(0); gt < STAR_DIINDEL::SIZE; ++gt)
        {
            dindel.phredLoghood[gt] = (unsigned) std::min(GermlineDiploidIndelSimpleGenotypeInfoCore::maxQ,
                                                          (int) (dindel.phredLoghood[gt] - minPL));
        }
    }

    return dindel;
}



/// convert the new AlleleGroupGenotype format to 0..N similar starling_diploid_indel intermediates as a
/// tempoary way for this method to communicate with the gVCF writer
///
/// TODO remove this function once we eliminate starling_diploid_indel as an intermediary format
static
bool
hackDiplotypeCallToCopyNumberCalls(
    const starling_base_options& opt,
    const starling_deriv_options& dopt,
    const pos_basecall_buffer& basecallBuffer,
    const pos_t targetPos,
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    const AlleleGroupGenotype& locusGenotype,
    const std::vector<LocusSupportingReadStats>& locusReadStats,
    const unsigned callerPloidy, ///< ploidy used by the genotyper
    const unsigned groupLocusPloidy, ///< actual ploidy in sample
    const bool isForcedOutput,
    gvcf_aggregator& gvcfer)
{
    static const bool is_tier2_pass(false);
    static const bool is_use_alt_indel(false);

    /// TODO STREL-125 generalize to multi-sample
    const unsigned sampleCount(1);

    bool isOutputAnyAlleles(false);

    // cycle through variant alleles in genotype:
    const unsigned alleleGroupSize(alleleGroup.size());
    for (unsigned genotypeAlleleIndex(0); genotypeAlleleIndex<alleleGroupSize; ++genotypeAlleleIndex)
    {
        if (not (isForcedOutput or (locusGenotype.maxGenotypeIndex != AG_GENOTYPE::HOMREF)))
        {
            continue;
        }

        const IndelKey& indelKey(alleleGroup.key(genotypeAlleleIndex));
        const IndelData& indelData(alleleGroup.data(genotypeAlleleIndex));

        if (indelKey.pos != targetPos) continue;

        assert ((not isForcedOutput) || indelData.isForcedOutput);

        isOutputAnyAlleles=true;

        starling_diploid_indel dindel(locusGenotypeToDindel(locusGenotype, genotypeAlleleIndex));

        // sample-independent info:
        std::unique_ptr<GermlineIndelLocusInfo> ii(new GermlineDiploidIndelLocusInfo(dopt.gvcf, sampleCount, indelKey, indelData, dindel));

        ii->anyVariantAlleleQuality = locusGenotype.anyVariantAlleleQuality;
        for (unsigned genotypeAlleleIndex2(0); genotypeAlleleIndex2<alleleGroupSize; ++genotypeAlleleIndex2)
        {
            const IndelKey& indelKey2(alleleGroup.key(genotypeAlleleIndex2));
            const IndelData& indelData2(alleleGroup.data(genotypeAlleleIndex2));
            ii->addAltIndelAllele(indelKey2, indelData2);
        }

        for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
        {
            auto& sampleInfo(ii->getSample(sampleIndex));
            sampleInfo.setPloidy(callerPloidy);
            if (callerPloidy != groupLocusPloidy)
            {
                sampleInfo.setPloidyConflict();
            }
            sampleInfo.supportCounts = locusReadStats[sampleIndex];

            // add info for PLs
            const unsigned fullAlleleCount(alleleGroupSize+1);
            auto& samplePLs(sampleInfo.genotypePhredLoghood);
            if (sampleInfo.getPloidy().isHaploid())
            {
                for (unsigned fullAlleleIndex(0); fullAlleleIndex<fullAlleleCount; ++fullAlleleIndex)
                {
                    samplePLs.getGenotypeLikelihood(fullAlleleIndex) = locusGenotype.phredLoghood[AG_GENOTYPE::getGenotypeId(fullAlleleIndex)];
                }
            }
            else if (sampleInfo.getPloidy().isDiploid())
            {
                for (unsigned fullAlleleIndex(0); fullAlleleIndex<fullAlleleCount; ++fullAlleleIndex)
                {
                    for (unsigned fullAlleleIndex2(fullAlleleIndex); fullAlleleIndex2<fullAlleleCount; ++fullAlleleIndex2)
                    {
                        samplePLs.getGenotypeLikelihood(fullAlleleIndex, fullAlleleIndex2) = locusGenotype.phredLoghood[AG_GENOTYPE::getGenotypeId(fullAlleleIndex, fullAlleleIndex2)];
                    }
                }
            }
            else
            {
                assert(false and "InvalidPloidy");
            }

            // add misc sample info from legacy sample indel report:
            const IndelSampleData& indelSampleData(indelData.getSampleData(sampleIndex));
            AlleleSampleReportInfo& indelSampleReportInfo(ii->getIndelSample(sampleIndex).reportInfo);
            getAlleleSampleReportInfo(opt, dopt, indelKey, indelSampleData, basecallBuffer,
                                      is_tier2_pass, is_use_alt_indel, indelSampleReportInfo);
        }

        gvcfer.add_indel(std::move(ii));
    }

    return isOutputAnyAlleles;
}



void
starling_pos_processor::
process_pos_indel_digt(const pos_t pos)
{
    /// TODO remove this legacy option, germline calling is *always and only* gvcf output now:
    assert(_opt.gvcf.is_gvcf_output());

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

    //   ************* end of sample generalization progress
    assert(sampleCount==1);
    const unsigned sampleIndex(0);
    sample_info& sif(sample(sampleIndex));

    // add candidate alleles to topVariant group, then order (including reference) according to
    // read evidence and take the top N allele candidates, N = callerPloidy
    //
    // track all forced output alleles in a separate group (even if they go into topVariant group)
    // to ensure that these are output even if not included in the most likely genotype
    //
    OrthogonalVariantAlleleCandidateGroup topVariantAlleleGroup;
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

        // rank and select top N alleles, N=callerPloidy
        assert(callerPloidy[sampleIndex]>0);
        selectTopOrthogonalAllelesInSample(sampleIndex, orthogonalVariantAlleles, callerPloidy[sampleIndex], topVariantAlleleGroup);
    }

    // at this point topVariantAlleleGroup represents the best alleles at the current position, now we
    // add conflicting alleles at other positions and re-rank, re-select the top alleles again:
    if (not topVariantAlleleGroup.empty())
    {
        addAllelesAtOtherPositions(pos, get_largest_total_indel_ref_span_per_read(), sampleIndex,
                                   callerPloidy[sampleIndex], getIndelBuffer(), topVariantAlleleGroup);
    }

    // genotype and report topVariantAlleleGroup
    //

    AlleleGroupGenotype locusGenotype;
    std::vector<LocusSupportingReadStats> locusReadStats(sampleCount);
    {
        // genotype the top N alleles
        getVariantAlleleGroupGenotypeLhoods(_opt, _dopt, sif.sample_opt, callerPloidy[sampleIndex], sampleIndex,
                                            topVariantAlleleGroup,
                                            locusGenotype, locusReadStats[sampleIndex]);

        // coerce output into older data-structures for gVCF output
        static const bool isForcedOutput(false);
        hackDiplotypeCallToCopyNumberCalls(_opt, _dopt, sif.bc_buff, pos, topVariantAlleleGroup,
                                           locusGenotype, locusReadStats, callerPloidy[sampleIndex],
                                           groupLocusPloidy[sampleIndex], isForcedOutput,
                                           *_gvcfer);
    }

    // score and report any remaining forced output alleles
    //
    {
        // trim the forced output allele set to take out any alleles already called as variants:
        if (not forcedOutputAlleleGroup.empty())
        {
            auto eraseForced = [&](const unsigned variantAlleleIndex)
            {
                const IndelKey& alleleKey(topVariantAlleleGroup.key(variantAlleleIndex));

                /// TMP: brute-force the gt match search:
                const unsigned forcedCount(forcedOutputAlleleGroup.size());
                for (unsigned forcedIndex(0); forcedIndex < forcedCount; ++forcedIndex)
                {
                    if (forcedOutputAlleleGroup.key(forcedIndex) == alleleKey)
                    {
                        forcedOutputAlleleGroup.alleles.erase(forcedOutputAlleleGroup.alleles.begin() + forcedIndex);
                        break;
                    }
                }
            };

            const unsigned topVariantAlleleCount(topVariantAlleleGroup.size());
            for (unsigned variantAlleleIndex(0); variantAlleleIndex < topVariantAlleleCount; variantAlleleIndex++)
            {
                if (AG_GENOTYPE::isAllelePresent(locusGenotype.maxGenotypeIndex, variantAlleleIndex))
                {
                    eraseForced(variantAlleleIndex);
                }
            }
        }

        // enumerate support for remaining forced output alleles compared to orthogonal genotyped variant alleles above
        const unsigned forcedOutputAlleleCount(forcedOutputAlleleGroup.size());
        for (unsigned forcedOutputAlleleIndex(0);
             forcedOutputAlleleIndex < forcedOutputAlleleCount; ++forcedOutputAlleleIndex)
        {
            AlleleGroupGenotype forcedAlleleLocusGenotype;
            getGenotypeLhoodsForForcedOutputAllele(_opt, _dopt, sif.sample_opt, callerPloidy[sampleIndex], sampleIndex,
                                                   topVariantAlleleGroup,
                                                   forcedOutputAlleleGroup, forcedOutputAlleleIndex,
                                                   forcedAlleleLocusGenotype, locusReadStats[sampleIndex]);

            // The above function should be set to force <*> and REF alleles into just REF for now,
            // so the most likely genotype should not contain the second allele
            //
            // This alleles compression should be relaxed once we have a way to express this in
            // the output.
            assert(not AG_GENOTYPE::isAllelePresent(forcedAlleleLocusGenotype.maxGenotypeIndex, 1));

            // fake an allele group with only the forced output allele so that we can output using
            // standard data structures
            OrthogonalVariantAlleleCandidateGroup fakeForcedOutputAlleleGroup;
            fakeForcedOutputAlleleGroup.addVariantAllele(forcedOutputAlleleGroup.alleles[forcedOutputAlleleIndex]);
            static const bool isForcedOutput(true);
            hackDiplotypeCallToCopyNumberCalls(_opt, _dopt, sif.bc_buff, pos, fakeForcedOutputAlleleGroup,
                                               forcedAlleleLocusGenotype, locusReadStats, callerPloidy[sampleIndex],
                                               groupLocusPloidy[sampleIndex], isForcedOutput, *_gvcfer);
        }
    }
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

        // sample-independent info:
        static const bool is_tier2_pass(false);
        static const bool is_use_alt_indel(true);

        bool isReportableAllele(isForcedOutput);

        std::unique_ptr<GermlineContinuousIndelLocusInfo> locusInfo(new GermlineContinuousIndelLocusInfo(sampleCount, pos));

        for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
        {
            locusInfo->getSample(sampleIndex).setPloidy(-1);

            const IndelSampleData& indelSampleData(indelData.getSampleData(sampleIndex));
            const sample_info& sif(sample(sampleIndex));

            AlleleSampleReportInfo& indelSampleReportInfo(locusInfo->getIndelSample(sampleIndex).reportInfo);
            getAlleleSampleReportInfo(_opt, _dopt, indelKey, indelSampleData, sif.bc_buff, is_tier2_pass,
                                      is_use_alt_indel, indelSampleReportInfo);

            if (not isReportableAllele)
            {
                if (indelSampleReportInfo.n_confident_indel_reads > 0) isReportableAllele = true;
            }
        }

        if (not isReportableAllele) continue;

        starling_continuous_variant_caller::add_indel_call(_opt, indelKey, indelData, *locusInfo);

        if (_opt.gvcf.is_gvcf_output())
        {
            _gvcfer->add_indel(std::move(locusInfo));
        }
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

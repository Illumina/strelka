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
#include "starling_common/OrthogonalVariantAlleleCandidateGroup.hh"
#include "starling_common/OrthogonalVariantAlleleCandidateGroupUtil.hh"
#include "starling_common/AlleleGroupGenotype.hh"
#include "blt_common/position_nonref_test.hh"
#include "blt_common/position_nonref_2allele_test.hh"
#include "blt_common/ref_context.hh"
#include "blt_util/log.hh"
#include "starling_continuous_variant_caller.hh"

#include <iomanip>
#include <starling_common/indel_util.hh>



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
    : base_t(opt,dopt,ref,streams,1),
      _opt(opt),
      _dopt(dopt),
      _streams(streams)
{
    static const unsigned sampleId(0);

    // setup gvcf aggregator
    if (_opt.gvcf.is_gvcf_output())
    {
        _gvcfer.reset(new gvcf_aggregator(
                          _opt,_dopt,ref,_nocompress_regions,
                          _streams.getSampleName(), _streams.gvcf_osptr(),
                          sample(sampleId).bc_buff));
    }

    // setup indel syncronizer:
    {
        sample_info& normal_sif(sample(0));

        double max_candidate_normal_sample_depth(-1.);
        if (dopt.gvcf.is_max_depth())
        {
            if (opt.max_candidate_indel_depth_factor > 0.)
            {
                max_candidate_normal_sample_depth = (opt.max_candidate_indel_depth_factor * dopt.gvcf.max_depth);
            }
        }

        if (opt.max_candidate_indel_depth > 0.)
        {
            if (max_candidate_normal_sample_depth > 0.)
            {
                max_candidate_normal_sample_depth = std::min(max_candidate_normal_sample_depth,static_cast<double>(opt.max_candidate_indel_depth));
            }
            else
            {
                max_candidate_normal_sample_depth = opt.max_candidate_indel_depth;
            }
        }

        const unsigned syncSampleId = getIndelBuffer().registerSample(normal_sif.estdepth_buff, normal_sif.estdepth_buff_tier2,
                                                                      max_candidate_normal_sample_depth);

        assert(syncSampleId == sampleId);

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
    if (sample_no!=0) return;

    sample_info& sif(sample(sample_no));

    const CleanedPileup& cpi(sif.cpi);
    const snp_pos_info& pi(cpi.rawPileup());

    _pileupCleaner.CleanPileupErrorProb(sif.cpi);

    const snp_pos_info& good_pi(cpi.cleanedPileup());
    const bool is_forced(is_forced_output_pos(pos));

    if (pi.calls.empty() && !is_forced) return;

    std::unique_ptr<GermlineSiteCallInfo> si(new GermlineContinuousSiteCallInfo(pos,pi.get_ref_base(),good_pi,
                                                                                _opt.used_allele_count_min_qscore, _opt.min_het_vf, is_forced));

    si->n_used_calls=cpi.n_used_calls();
    si->n_unused_calls=cpi.n_unused_calls();
    // hpol filter
    si->hpol=get_snp_hpol_size(pos,_ref);


    starling_continuous_variant_caller::position_snp_call_continuous(_opt, good_pi, (GermlineContinuousSiteCallInfo&)*si);

    if (_opt.is_counts)
    {
        report_counts(good_pi,si->n_unused_calls,si->pos+1,*_streams.counts_osptr());
    }
    if (si->is_snp())
    {
        // this needs to be updated no matter where the snp-call is written to:
        if (_is_variant_windows)
        {
            _variant_print_pos.insert(pos);
            _is_skip_process_pos=false;
        }
    }

    if (_opt.gvcf.is_gvcf_output())
    {
        _gvcfer->add_site(std::move(si));
    }
}



void
starling_pos_processor::
process_pos_snp_single_sample_impl(
    const pos_t pos,
    const unsigned sample_no)
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

    sample_info& sif(sample(sample_no));

    const CleanedPileup& cpi(sif.cpi);
    const snp_pos_info& pi(cpi.rawPileup());

    // note multi-sample status -- can still be called only for one sample
    // and only for sample 0. working on generalization:
    //
    if (sample_no!=0) return;

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


    std::unique_ptr<GermlineDiploidSiteCallInfo> si(new GermlineDiploidSiteCallInfo(pos,pi.get_ref_base(),good_pi,_opt.used_allele_count_min_qscore, is_forced));
    si->n_used_calls=cpi.n_used_calls();
    si->n_unused_calls=cpi.n_unused_calls();


    // delay writing any snpcalls so that anomaly tests can (optionally) be applied as filters:
    //
    nonref_test_call nrc;
    //lrt_snp_call lsc;
    //monoploid_genotype mgt;
    //std::unique_ptr<nploid_genotype> ngt_ptr;

    // check whether we're in a haploid region:
    si->dgt.ploidy=(get_ploidy(pos));

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
    if (_opt.is_bsnp_monoploid)
    {
        position_snp_call_pprob_monogt(_opt.bsnp_monoploid_theta,good_pi,mgt);
    }
    if (_opt.is_bsnp_nploid)
    {
        ngt_ptr.reset(new nploid_genotype(*_ninfo));
        position_snp_call_pprob_nploid(_opt.bsnp_nploid_snp_prob,good_pi,*_ninfo,*ngt_ptr);
    }
#endif

    //    const bool is_snp(nrc.is_snp || lsc.is_snp || _site_info.dgt.is_snp || mgt.is_snp || (ngt_ptr.get() && ngt_ptr->is_snp));
    const bool is_snp(nrc.is_snp || si->dgt.is_snp);

    // find anomalies:
    //
#if 0
    bool is_pos_adis(false);
    bool is_pos_acov(false);

    if ((_opt.is_adis_table || _opt.is_adis_lrt) && is_snp)
    {
        if (_opt.is_adis_table)
        {
            is_pos_adis = (is_pos_adis || position_strand_distro_anomaly(_opt.adis_table_alpha,good_pi,_ws));
        }
        if (_opt.is_adis_lrt)
        {
            is_pos_adis = (is_pos_adis || position_strand_distro_anomaly_lrt(_opt.adis_lrt_alpha,good_pi));
        }
    }
    if (_opt.is_acov)
    {
        is_pos_acov = position_strand_coverage_anomaly(_opt.acov_alpha,pi);
    }
#endif

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
#endif
        if (si->dgt.is_snp)
        {
            // this needs to be updated no matter where the snp-call is written to:
            if (_is_variant_windows)
            {
                _variant_print_pos.insert(pos);
                _is_skip_process_pos=false;
            }
        }
#if 0
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
process_pos_indel_single_sample(
    const pos_t pos,
    const unsigned sample_no)
{
    if (_opt.is_bsnp_diploid())
    {
        process_pos_indel_single_sample_digt(pos,sample_no);
    }
    else
    {
        process_pos_indel_single_sample_continuous(pos,sample_no);
    }
}



/// TEMPORARY
/// refine topVariantAlleleGroup to consider alts shared by all top alleles, then rerank
/// and re-select top groupLocusPloidy alleles
static
void
addAllelesAtOtherPositions(
    const pos_t pos,
    const pos_t largest_total_indel_ref_span_per_read,
    const unsigned sampleId,
    const int groupLocusPloidy,
    const IndelBuffer& indelBuffer,
    OrthogonalVariantAlleleCandidateGroup& alleleGroup)
{
    const pos_t minIndelBufferPos(pos-largest_total_indel_ref_span_per_read);

    const unsigned inputAlleleCount(alleleGroup.size());

    // first get the set of candidate alt alleles from another position
    //
    std::vector<IndelKey> filteredAltAlleles;
    {
        const known_pos_range inputAlleleGroupRange(alleleGroup.getReferenceRange());

        // extend end_pos by one to ensure that we find indels adjacent to the right end of the range
        /// TODO add strict definition and unit tests to rangeIterator wrt adjacent indels
        const auto indelIterPair(indelBuffer.rangeIterator(inputAlleleGroupRange.begin_pos, inputAlleleGroupRange.end_pos+1));
        for (auto altAlleleIter(indelIterPair.first); altAlleleIter!=indelIterPair.second; ++altAlleleIter)
        {
            const IndelKey& altAlleleKey(altAlleleIter->first);

            // all alleles with this starting position have already been
            // considered..
            if (altAlleleKey.pos == pos) continue;

            // filter out indels which have already been cleared out of the indel buffer:
            if (altAlleleKey.pos < minIndelBufferPos) continue;

            // no breakpoints:
            if (altAlleleKey.is_breakpoint()) continue;

            // must be orthogonal to all input alleles:
            {
                bool isOrthogonalToAllInputAlleles(true);
                for (unsigned inputAlleleIndex(0); inputAlleleIndex < inputAlleleCount; inputAlleleIndex++)
                {
                    const IndelKey& inputAlleleKey(alleleGroup.key(inputAlleleIndex));
                    if (is_indel_conflict(altAlleleKey, inputAlleleKey)) continue;
                    isOrthogonalToAllInputAlleles = false;
                    break;
                }
                if (not isOrthogonalToAllInputAlleles) continue;
            }

            const IndelData& altAlleleData(getIndelData(altAlleleIter));

            // must be a candidate allele:
            if (not indelBuffer.isCandidateIndel(altAlleleKey, altAlleleData)) continue;

            // made it! add allele to the set we move forward with:
            filteredAltAlleles.push_back(altAlleleKey);
        }
    }

    if (filteredAltAlleles.empty()) return;

    OrthogonalVariantAlleleCandidateGroup altAlleleGroup;
    for (const auto& altAlleleKey : filteredAltAlleles)
    {
        const auto altIter(indelBuffer.getIndelIter(altAlleleKey));
        altAlleleGroup.addVariantAllele(altIter);
    }

    if (altAlleleGroup.size()>1)
    {
        // rank alt alleles and include from highest to lowest unless interference clique is broken:
        unsigned referenceRank(0);
        OrthogonalVariantAlleleCandidateGroup rankedAltAlleleGroup;
        rankOrthogonalAllelesInSample(sampleId, altAlleleGroup, rankedAltAlleleGroup, referenceRank);

        altAlleleGroup.clear();

        for (const auto& rankedAltAlleleIter : rankedAltAlleleGroup.alleles)
        {
            bool isGroupOrthogonal(true);
            for (const auto& altAlleleIter : altAlleleGroup.alleles)
            {
                if (not is_indel_conflict(rankedAltAlleleIter->first, altAlleleIter->first))
                {
                    isGroupOrthogonal = false;
                    break;
                }
            }
            if (isGroupOrthogonal)
            {
                altAlleleGroup.alleles.push_back(rankedAltAlleleIter);
            }
        }
    }

    // put all qualifying alts back together with variants to form an extended allele set:
    OrthogonalVariantAlleleCandidateGroup extendedVariantAlleleGroup(alleleGroup);
    for (const auto& altAlleleIter : altAlleleGroup.alleles)
    {
        extendedVariantAlleleGroup.alleles.push_back(altAlleleIter);
    }

    // rerank and reselect top N alleles, N=groupLocusPloidy
    //
    selectTopOrthogonalAllelesInSample(sampleId, extendedVariantAlleleGroup, groupLocusPloidy, alleleGroup);
}



/// TEMPORARY
static
void
insertIndelInGvcf(
    const starling_base_options& opt,
    const starling_base_deriv_options& dopt,
    const reference_contig_segment& ref,
    const pos_basecall_buffer& basecallBuffer,
    const IndelKey& indelKey,
    const IndelData& indelData,
    const unsigned sampleId,
    const starling_diploid_indel& dindel,
    gvcf_aggregator& gvcfer)
{
    static const bool is_tier2_pass(false);
    static const bool is_use_alt_indel(false);

    // sample-independent info:
    starling_indel_report_info indelReportInfo;
    get_starling_indel_report_info(indelKey, indelData, ref, indelReportInfo);

    // sample-specific info: (division doesn't really matter
    // in single-sample case)
    const IndelSampleData& indelSampleData(indelData.getSampleData(sampleId));
    starling_indel_sample_report_info indelSampleReportInfo;
    get_starling_indel_sample_report_info(opt, dopt, indelKey, indelSampleData, basecallBuffer,
                                          is_tier2_pass, is_use_alt_indel, indelSampleReportInfo);

    gvcfer.add_indel(std::unique_ptr<GermlineIndelCallInfo>(new GermlineDiploidIndelCallInfo(indelKey, indelData, dindel, indelReportInfo, indelSampleReportInfo)));
}



/// TEMPORARY
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
            switch(static_cast<index_t>(genotypeId))
            {
                case HOM0: return STAR_DIINDEL::HOM;
                case HET0:
                case HET01: return STAR_DIINDEL::HET;
                default:    return STAR_DIINDEL::NOINDEL;
            }
        }
        else if (alleleId == 1)
        {
            switch(static_cast<index_t>(genotypeId))
            {
                case HOM1: return STAR_DIINDEL::HOM;
                case HET1:
                case HET01: return STAR_DIINDEL::HET;
                default:    return STAR_DIINDEL::NOINDEL;
            }
        }
        else
        {
            assert(false and "Unknown alleleId");
            return STAR_DIINDEL::SIZE;
        }
    }
}



/// TEMPORARY
static
starling_diploid_indel
locusGenotypeToDindel(
    const AlleleGroupGenotype& locusGenotype,
    const unsigned alleleId,
    const unsigned groupLocusPloidy,
    const bool isForcedOuput,
    const bool isZeroCoverage)
{
    starling_diploid_indel dindel;
    dindel.is_forced_output = isForcedOuput;
    dindel.is_zero_coverage = isZeroCoverage;
    dindel.ploidy = groupLocusPloidy;

    dindel.pprob[STAR_DIINDEL::NOINDEL] = locusGenotype.posteriorProb[AG_GENOTYPE::HOMREF];
    dindel.pprob[STAR_DIINDEL::HET] = locusGenotype.posteriorProb[AG_GENOTYPE::getAlleleHetId(alleleId)];
    dindel.pprob[STAR_DIINDEL::HOM] = locusGenotype.posteriorProb[AG_GENOTYPE::getAlleleHomId(alleleId)];

    if (locusGenotype.maxGenotypeIndex == AG_GENOTYPE::HET01)
    {
        dindel.pprob[STAR_DIINDEL::HET] = locusGenotype.posteriorProb[AG_GENOTYPE::HET01];
    }

    dindel.indel_qphred = (int) locusGenotype.variantAlleleQuality;
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



/// TEMPORARY
///
/// Note that ploidy provided here is the actual ploidy value for this locus rather than the one used by the
/// caller, these two values should only be different when ploidy==0
///
static
bool
hackDiplotypeCallToCopyNumberCalls(
    const starling_base_options& opt,
    const starling_base_deriv_options& dopt,
    const reference_contig_segment& ref,
    const pos_basecall_buffer& basecallBuffer,
    const unsigned sampleId,
    const pos_t targetPos,
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    const AlleleGroupGenotype& locusGenotype,
    const unsigned groupLocusPloidy,
    const bool isForcedOutput,
    gvcf_aggregator& gvcfer)
{
    bool isOutputAnyAlleles(false);

    // cycle through variant alleles in genotype:
    const unsigned alleleGroupSize(alleleGroup.size());
    for (unsigned genotypeAlleleIndex(0); genotypeAlleleIndex<alleleGroupSize; ++genotypeAlleleIndex)
    {
        if (not (isForcedOutput or
                 AG_GENOTYPE::isAllelePresent(locusGenotype.maxGenotypeIndex, genotypeAlleleIndex)))
        {
            continue;
        }

        const IndelKey& indelKey(alleleGroup.key(genotypeAlleleIndex));
        const IndelData& indelData(alleleGroup.data(genotypeAlleleIndex));

        if (indelKey.pos != targetPos) continue;

        isOutputAnyAlleles=true;

        const IndelSampleData& indelSampleData(indelData.getSampleData(sampleId));
        const bool isZeroCoverage(indelSampleData.read_path_lnp.empty());

        starling_diploid_indel dindel(locusGenotypeToDindel(locusGenotype, genotypeAlleleIndex,
                                                            groupLocusPloidy, isForcedOutput, isZeroCoverage));
        insertIndelInGvcf(opt, dopt, ref, basecallBuffer, indelKey, indelData, sampleId, dindel, gvcfer);
    }
    return isOutputAnyAlleles;
}



void
starling_pos_processor::
process_pos_indel_single_sample_digt(
    const pos_t pos,
    const unsigned sampleId)
{
    // note multi-sample status -- can still be called only for one sample
    // and only for sample 0. working on generalization:
    //
    if (sampleId!=0) return;

    /// TODO remove this legacy option, germline calling is *always and only* gvcf output now:
    assert(_opt.gvcf.is_gvcf_output());

    sample_info& sif(sample(sampleId));
    auto it(getIndelBuffer().positionIterator(pos));
    const auto it_end(getIndelBuffer().positionIterator(pos + 1));

    // define groups of overlapping alleles:
    //
    // alleles should form "conflict graphs", where an edge exists
    // between two alleles which cannot exist together on the same haplotype
    //
    // When only evaluating alleles at a single position in one sample, this graph is a clique,
    // so the process of picking the mostly likely candidates is simple, we can pick the top
    // two alleles (by some evidence metric), and add in the refernece if it is not one of the
    // top two.
    //
    // This method will have to be redesigned for multi-sample and/or non-clique conflict graphs.
    //
    OrthogonalVariantAlleleCandidateGroup orthogonalVariantAlleles;
    for (; it!=it_end; ++it)
    {
        const IndelKey& indelKey(it->first);
        const IndelData& indelData(getIndelData(it));

        if (indelKey.is_breakpoint()) continue;

        const bool isForcedOutput(indelData.is_forced_output);
        if (not isForcedOutput)
        {
            const IndelSampleData& indelSampleData(indelData.getSampleData(sampleId));
            const bool isZeroCoverage(indelSampleData.read_path_lnp.empty());

            if (isZeroCoverage) continue;
            if (not getIndelBuffer().isCandidateIndel(indelKey, indelData)) continue;
        }

        // all alleles at the same position are automatically conflicting/orthogonal:
        orthogonalVariantAlleles.addVariantAllele(it);
    }

    if (orthogonalVariantAlleles.empty()) return;

    // determine ploidy for this locus
    //
    // Assume entire allele group is covered by one ploidy type in nearly all cases,
    // in case of a conflict use the highest ploidy overlapped by the group.
    //
    unsigned groupLocusPloidy(0);
    {
        const known_pos_range alleleGroupRange(orthogonalVariantAlleles.getReferenceRange());
        const unsigned groupLeftPloidy(get_ploidy(alleleGroupRange.begin_pos));
        const unsigned groupRightPloidy(get_ploidy(alleleGroupRange.end_pos));

        groupLocusPloidy=std::max(groupLeftPloidy,groupRightPloidy);
    }

    // groupLocusPloidy of 0 is treated as a special case, if this happens the
    // entire calling method reverts to a ploidy of 2, but the locus ploidy is
    // passed into the gVCF writer as 0. The gVCF writer can decide what to do
    // with this information from there.
    //
    const unsigned callerPloidy((groupLocusPloidy==0) ? 2 : groupLocusPloidy);

    // if callerPloidy > orthogonalVariantAlleleCount, order alleles according to read evidence,
    // and take the top 'callerPloidy' allele candidates (besides the reference)
    //
    // track all forcedOutput alleles in a separate group (even if they go into topVariant group)
    // to ensure that these are output even if not included in the most likely genotype
    //
    OrthogonalVariantAlleleCandidateGroup topVariantAlleleGroup;
    OrthogonalVariantAlleleCandidateGroup forcedOutputAlleleGroup;
    {
        const unsigned orthogonalVariantAlleleCount(orthogonalVariantAlleles.size());
        for (unsigned alleleIndex(0); alleleIndex < orthogonalVariantAlleleCount; alleleIndex++)
        {
            const IndelData& indelData(orthogonalVariantAlleles.data(alleleIndex));
            if (indelData.is_forced_output)
            {
                forcedOutputAlleleGroup.addVariantAllele(orthogonalVariantAlleles.iter(alleleIndex));
            }
        }

        // rank and select top N alleles, N=callerPloidy
        assert(callerPloidy>0);
        selectTopOrthogonalAllelesInSample(sampleId, orthogonalVariantAlleles, callerPloidy, topVariantAlleleGroup);
    }

    // refine topVariantAlleleGroup to consider alts shared by all top alleles, then rerank
    if (not topVariantAlleleGroup.empty())
    {
        addAllelesAtOtherPositions(pos, get_largest_total_indel_ref_span_per_read(), sampleId,
                                   callerPloidy, getIndelBuffer(), topVariantAlleleGroup);
    }

    // score and report topVariantAlleleGroup
    //
    AlleleGroupGenotype locusGenotype;
    bool isOutputAlleles(false);
    {
        getVariantAlleleGroupGenotypeLhoods(_dopt, sif.sample_opt, _dopt.getIndelGenotypePriors(), callerPloidy,
                                            sampleId, topVariantAlleleGroup, locusGenotype);

        // (1) TEMPORARY hack called variant alleles into old dindel structure(s)
        static const bool isForcedOutput(false);
        isOutputAlleles = (hackDiplotypeCallToCopyNumberCalls(_opt, _dopt, _ref, sif.bc_buff, sampleId, pos,
                                                              topVariantAlleleGroup, locusGenotype,
                                                              groupLocusPloidy, isForcedOutput, *_gvcfer));
    }

    // score and report any remaining forced output alleles
    //
    {
        /// (2) trim the forcedGT allele set to take out any alleles already called as variants:
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

        if (not forcedOutputAlleleGroup.empty())
        {
            const unsigned topVariantAlleleCount(topVariantAlleleGroup.size());
            for (unsigned variantAlleleIndex(0); variantAlleleIndex < topVariantAlleleCount; variantAlleleIndex++)
            {
                if (AG_GENOTYPE::isAllelePresent(locusGenotype.maxGenotypeIndex, variantAlleleIndex))
                {
                    eraseForced(variantAlleleIndex);
                }
            }
        }

        // enumerate support for remaining forcedOutput alleles compared to orthogonal genotyped variant alleles above
        const unsigned forcedOutputAlleleCount(forcedOutputAlleleGroup.size());
        for (unsigned forcedOutputAlleleIndex(0);
             forcedOutputAlleleIndex < forcedOutputAlleleCount; ++forcedOutputAlleleIndex)
        {
            AlleleGroupGenotype forcedAlleleLocusGenotype;
            getGenotypeLhoodsForForcedOutputAllele(_dopt, sif.sample_opt, _dopt.getIndelGenotypePriors(),
                                                   callerPloidy, sampleId, topVariantAlleleGroup,
                                                   forcedOutputAlleleGroup,
                                                   forcedOutputAlleleIndex, forcedAlleleLocusGenotype);

            // TEMPORARY, the above function should be set to compress <*> and REF alleles to just REF for now,
            // so the most likely genotype should not contain the second allele
            assert(not AG_GENOTYPE::isAllelePresent(forcedAlleleLocusGenotype.maxGenotypeIndex, 1));

            // TEMPORARY fake an allele group with only the forced output allele so that we can output using
            // standard data structures
            OrthogonalVariantAlleleCandidateGroup fakeForcedOutputAlleleGroup;
            fakeForcedOutputAlleleGroup.addVariantAllele(forcedOutputAlleleGroup.alleles[forcedOutputAlleleIndex]);
            static const bool isForcedOutput(true);
            hackDiplotypeCallToCopyNumberCalls(_opt, _dopt, _ref, sif.bc_buff, sampleId,
                                               pos, fakeForcedOutputAlleleGroup,
                                               forcedAlleleLocusGenotype, groupLocusPloidy, isForcedOutput,
                                               *_gvcfer);

            isOutputAlleles = true;
        }
    }

    // finalize response to any reported indels for this pos:
    if (isOutputAlleles and _is_variant_windows)
    {
        _variant_print_pos.insert(pos);
        _is_skip_process_pos=false;
    }

#if 0
    // TODO implement indel overlap resolution
    //
    // punt conflict resolution for now....
    {
        // indel_report_info needs to be run first now so that
        // local small repeat info is available to the indel
        // caller

        // sample-independent info:
        starling_indel_report_info indelReportInfo;
        get_starling_indel_report_info(indelKey,indelData,_ref,indelReportInfo);

        // STARKA-248 filter invalid indel
        /// TODO: filter this issue earlier (occurs as, e.g. 1D1I which matches ref)
        if (indelReportInfo.vcf_indel_seq == indelReportInfo.vcf_ref_seq) continue;

        static const bool is_tier2_pass(false);
        static const bool is_use_alt_indel(true);

        starling_diploid_indel dindel;
        dindel.is_forced_output = isForcedOutput;
        dindel.is_zero_coverage = isZeroCoverage;

        {
            // check whether we're in a haploid/noploid region, for indels just check
            // start position and end position, approximating that the whole
            // region in between has the same ploidy, for any anomalous state
            // revert to 'noploid':
            const int indelLeftPloidy(get_ploidy(indelKey.pos));
            const int indelRightPloidy(get_ploidy(indelKey.right_pos()));

            if (indelLeftPloidy == indelRightPloidy)
            {
                dindel.ploidy = indelLeftPloidy;
            }
            else
            {
                dindel.ploidy = 0;
            }
        }

        _dopt.incaller().starling_indel_call_pprob_digt(
            _opt,_dopt,
            sif.sample_opt,
            indelKey,indelSampleData,is_use_alt_indel,dindel);

        bool is_indel(false);
        if ((dindel.is_indel) || (dindel.is_forced_output))
        {
            is_indel=true;

            // sample-specific info: (division doesn't really matter
            // in single-sample case)
            starling_indel_sample_report_info isri;
            get_starling_indel_sample_report_info(_opt, _dopt,indelKey,indelSampleData,sif.bc_buff,
                                                  is_tier2_pass,is_use_alt_indel,isri);

            if (_opt.gvcf.is_gvcf_output())
            {
                assert(indelKey.pos==pos);
                _gvcfer->add_indel(std::unique_ptr<GermlineIndelCallInfo>(new GermlineDiploidIndelCallInfo(indelKey,indelData, dindel,indelReportInfo,isri)));
            }

            if (_is_variant_windows)
            {
                _variant_print_pos.insert(pos);
                _is_skip_process_pos=false;
            }
        }

        /// \TODO put this option under runtime control...
        /// \TODO setup option so that read keys persist longer when needed for this case...
        ///
        static const bool is_print_indel_evidence(false);

        if (is_print_indel_evidence && is_indel)
        {
            std::ostream& report_os(std::cerr);
            report_os << "INDEL_EVIDENCE " << indelKey;

            for (const auto& val : indelSampleData.read_path_lnp)
            {
                const align_id_t read_id(val.first);
                const ReadPathScores& lnp(val.second);
                const ReadPathScores pprob(indel_lnp_to_pprob(_dopt,lnp,is_tier2_pass,is_use_alt_indel));
                const starling_read* srptr(sif.read_buff.get_read(read_id));

                report_os << "read key: ";
                if (NULL==srptr) report_os << "UNKNOWN_KEY";
                else            report_os << srptr->key();
                report_os << "\n"
                          << "read log_lhoods: " << lnp << "\n"
                          << "read pprobs: " << pprob << "\n";
            }
        }
    }
#endif
}



void
starling_pos_processor::
process_pos_indel_single_sample_continuous(
    const pos_t pos,
    const unsigned sampleId)
{
    // note multi-sample status -- can still be called only for one sample
    // and only for sample 0. working on generalization:
    //
    if (sampleId!=0) return;

    sample_info& sif(sample(sampleId));
    auto it(getIndelBuffer().positionIterator(pos));
    const auto it_end(getIndelBuffer().positionIterator(pos + 1));

    std::unique_ptr<GermlineContinuousIndelCallInfo> info;

    for (; it!=it_end; ++it)
    {
        const IndelKey& indelKey(it->first);
        const IndelData& indelData(getIndelData(it));
        const bool isForcedOutput(indelData.is_forced_output);

        const IndelSampleData& indelSampleData(indelData.getSampleData(sampleId));
        const bool isZeroCoverage(indelSampleData.read_path_lnp.empty());

        if (! isForcedOutput)
        {
            if (isZeroCoverage) continue;
            if (!getIndelBuffer().isCandidateIndel(indelKey, indelData)) continue;
        }

        // sample-independent info:
        starling_indel_report_info indelReportInfo;
        get_starling_indel_report_info(indelKey,indelData,_ref,indelReportInfo);

        static const bool is_tier2_pass(false);
        static const bool is_use_alt_indel(true);

        if (!info)
            info.reset(new GermlineContinuousIndelCallInfo(pos));

        starling_indel_sample_report_info isri;
        get_starling_indel_sample_report_info(_opt, _dopt,indelKey,indelSampleData,sif.bc_buff, is_tier2_pass,is_use_alt_indel,isri);
        starling_continuous_variant_caller::add_indel_call(_opt, indelKey, indelData, indelReportInfo, isri, *info);
    }
    if (info && (info->is_indel() || info->is_forced_output()))
    {
        if (_opt.gvcf.is_gvcf_output())
        {
            _gvcfer->add_indel(std::move(info));
        }

        if (_is_variant_windows)
        {
            _variant_print_pos.insert(pos);
            _is_skip_process_pos=false;
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

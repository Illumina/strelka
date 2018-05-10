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

#include "position_somatic_snv_strand_grid.hh"
#include "position_somatic_snv_strand_grid_vcf.hh"
#include "somatic_indel_grid.hh"
#include "strelka_pos_processor.hh"
#include "blt_util/log.hh"
#include "starling_common/AlleleReportInfoUtil.hh"
#include "starling_common/starling_pos_processor_base_stages.hh"

#include <iomanip>



strelka_pos_processor::
strelka_pos_processor(
    const strelka_options& opt,
    const strelka_deriv_options& dopt,
    const reference_contig_segment& ref,
    const strelka_streams& fileStreams,
    RunStatsManager& statsManager)
    : base_t(opt, dopt, ref, fileStreams, STRELKA_SAMPLE_TYPE::SIZE, statsManager)
    , _opt(opt)
    , _dopt(dopt)
    , _streams(fileStreams)
    , _scallProcessor(fileStreams.somatic_callable_osptr())
    , _indelWriter(opt, dopt, fileStreams.somatic_indel_osptr())
    , _indelRegionIndexNormal(0)
    , _indelRegionIndexTumor(0)
{
    using namespace STRELKA_SAMPLE_TYPE;

    sample_info& normal_sif(sample(NORMAL));
    sample_info& tumor_sif(sample(TUMOR));

    // set sample-specific parameter overrides:
    normal_sif.sampleOptions.min_read_bp_flank = opt.normal_sample_min_read_bp_flank;

    // setup indel buffer samples:
    {
        /// TODO: setup a stronger sample id handler -- using the version from manta would be a good start here
        sample_id_t sample_id;
        sample_id = getIndelBuffer().registerSample(normal_sif.estdepth_buff, normal_sif.estdepth_buff_tier2, true);

        assert(sample_id == NORMAL);

        sample_id = getIndelBuffer().registerSample(tumor_sif.estdepth_buff, tumor_sif.estdepth_buff_tier2, false);

        assert(sample_id == TUMOR);

        getIndelBuffer().finalizeSamples();
    }

    // setup indel avg window:
    _indelRegionIndexNormal= normal_sif.localRegionStatsCollection.addNewLocalRegionStatsSize(opt.sfilter.indelRegionFlankSize * 2);
    _indelRegionIndexTumor= tumor_sif.localRegionStatsCollection.addNewLocalRegionStatsSize(opt.sfilter.indelRegionFlankSize * 2);
}



void
strelka_pos_processor::
reset()
{
    base_t::reset();

    _indelWriter.clear();
    _noisePos.clear();
}



void
strelka_pos_processor::
resetRegion(
    const std::string& chromName,
    const known_pos_range2& regionRange)
{
    base_t::resetRegionBase(chromName, regionRange);

    // setup norm and max filtration depths
    {
        if (_dopt.sfilter.is_max_depth())
        {
            /// TODO verify that chroms match bam chroms
            cdmap_t::const_iterator cdi(_dopt.sfilter.chrom_depth.find(_chromName));
            if (cdi == _dopt.sfilter.chrom_depth.end())
            {
                std::ostringstream oss;
                oss << "Can't find chromosome: '" << _chromName << "' in chrom depth file: '"
                    << _opt.sfilter.chrom_depth_file << "'";
                throw blt_exception(oss.str().c_str());
            }
            _normChromDepth = cdi->second;
            _maxChromDepth = (_normChromDepth * _opt.sfilter.max_depth_factor);
        }
        assert(_normChromDepth >= 0.);
        assert(_maxChromDepth >= 0.);
    }

    // setup indel buffer max depth:
    {
        double maxIndelCandidateDepthSumOverNormalSamples(-1.);
        if (_dopt.sfilter.is_max_depth())
        {
            if (_opt.max_candidate_indel_depth_factor > 0.)
            {
                maxIndelCandidateDepthSumOverNormalSamples = (_opt.max_candidate_indel_depth_factor * _maxChromDepth);
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
strelka_pos_processor::
insert_noise_pos(
    const pos_t pos,
    const SiteNoise& sn)
{
    _stagemanPtr->validate_new_pos_value(pos,STAGE::READ_BUFFER);
    _noisePos.insertSiteNoise(pos,sn);
}



void
strelka_pos_processor::
process_pos_snp_somatic(const pos_t pos)
{
    using namespace STRELKA_SAMPLE_TYPE;

    const pos_t output_pos(pos+1);

    sample_info& normal_sif(sample(NORMAL));
    sample_info& tumor_sif(sample(TUMOR));

    // TODO this is ridiculous -- if the tier2 data scheme works then come back and clean this up:
    static const unsigned n_tier(2);
    CleanedPileup* normal_cpi_ptr[n_tier] = { &(normal_sif.cleanedPileup), &(_tier2_cpi[NORMAL]) };
    CleanedPileup* tumor_cpi_ptr[n_tier] = { &(tumor_sif.cleanedPileup), &(_tier2_cpi[TUMOR]) };

    for (unsigned t(0); t<n_tier; ++t)
    {
        const bool is_include_tier2(t!=0);
        if (is_include_tier2 && (! _opt.useTier2Evidence)) continue;
        _pileupCleaner.CleanPileup(normal_sif.basecallBuffer.get_pos(pos),is_include_tier2,*(normal_cpi_ptr[t]));
        _pileupCleaner.CleanPileup(tumor_sif.basecallBuffer.get_pos(pos),is_include_tier2,*(tumor_cpi_ptr[t]));
    }

    // note single-sample anomaly filtration won't apply here (more of
    // a vestigial blt feature anyway)
    //

    // retain original blt loop structure from the single-sample case
    // to allow for multiple interacting tests at one site
    //

    //    somatic_snv_genotype sgt;
    somatic_snv_genotype_grid sgtg;

    if (_opt.is_somatic_snv())
    {
        sgtg.is_forced_output=is_forced_output_pos(pos);

        const extended_pos_info* normal_epi_t2_ptr(nullptr);
        const extended_pos_info* tumor_epi_t2_ptr(nullptr);
        if (_opt.useTier2Evidence)
        {
            normal_epi_t2_ptr=&(normal_cpi_ptr[1]->getExtendedPosInfo());
            tumor_epi_t2_ptr=&(tumor_cpi_ptr[1]->getExtendedPosInfo());
        }

        const bool isComputeNonSomatic(_opt.is_somatic_callable());

        _dopt.sscaller_strand_grid().position_somatic_snv_call(
            normal_cpi_ptr[0]->getExtendedPosInfo(),
            tumor_cpi_ptr[0]->getExtendedPosInfo(),
            normal_epi_t2_ptr,
            tumor_epi_t2_ptr,
            isComputeNonSomatic,
            sgtg);

        if (_opt.is_somatic_callable())
        {
            _scallProcessor.addToRegion(_chromName,output_pos,sgtg);
        }
    }

    // report events:
    //
    if (sgtg.is_output())
    {
        {
            const SiteNoise* snp(_noisePos.getPos(pos));
            if (snp == nullptr)
            {
                sgtg.sn.clear();
            }
            else
            {
                sgtg.sn = *snp;
            }
        }
        std::ostream& bos(*_streams.somatic_snv_osptr());

        // have to keep tier1 counts for filtration purposes:
#ifdef SOMATIC_DEBUG
        write_snv_prefix_info_file(_chromName,output_pos,ref_base,normald,tumord,log_os);
        log_os << "\n";
#endif

        bos << _chromName << '\t'
            << output_pos << '\t'
            << ".";

        static const bool is_write_nqss(false);
        write_vcf_somatic_snv_genotype_strand_grid(_opt, _dopt, sgtg, is_write_nqss, *(normal_cpi_ptr[0]),
                                                   *(tumor_cpi_ptr[0]), *(normal_cpi_ptr[1]), *(tumor_cpi_ptr[1]),
                                                   _normChromDepth, _maxChromDepth, bos);
        bos << "\n";
    }
}



void
strelka_pos_processor::
process_pos_variants_impl(
    const pos_t pos,
    const bool isPosPrecedingReportableRange)
{
    if (isPosPrecedingReportableRange) return;

    try
    {
        process_pos_snp_somatic(pos);
    }
    catch (...)
    {
        log_os << "Exception caught while attempting to call somatic SNV at position: " << (pos+1) << "\n";
        throw;
    }

    try
    {
        process_pos_indel_somatic(pos);
    }
    catch (...)
    {
        log_os << "Exception caught while attempting to call somatic indel at position: " << (pos+1) << "\n";
        throw;
    }
}



void
strelka_pos_processor::
process_pos_indel_somatic(const pos_t pos)
{
    using namespace STRELKA_SAMPLE_TYPE;

    //    std::ostream& report_os(get_report_os());
    sample_info& normal_sif(sample(NORMAL));
    sample_info& tumor_sif(sample(TUMOR));

    auto indelIter(getIndelBuffer().positionIterator(pos));
    const auto indelIterEnd(getIndelBuffer().positionIterator(pos + 1));

    for (; indelIter != indelIterEnd; ++indelIter)
    {
        const IndelKey& indelKey(indelIter->first);

        // don't write breakpoint output:
        if (indelKey.is_breakpoint()) continue;

        if (indelKey.isMismatch()) continue;

        const IndelData& indelData(getIndelData(indelIter));

        if (!getIndelBuffer().isCandidateIndel(indelKey, indelData)) continue;

        const IndelSampleData& normalIndelSampleData(indelData.getSampleData(NORMAL));
        const IndelSampleData& tumorIndelSampleData(indelData.getSampleData(TUMOR));

        if (not indelData.isForcedOutput)
        {
            if (normalIndelSampleData.read_path_lnp.empty() && tumorIndelSampleData.read_path_lnp.empty()) continue;
        }

        if (_opt.is_somatic_indel())
        {
            // indel_report_info needs to be run first now so that
            // local small repeat info is available to the indel
            // caller
            // indel summary info
            SomaticIndelVcfInfo siInfo;

            getSingleIndelAlleleVcfSummaryStrings(indelKey, indelData, _ref, siInfo.vcf_indel_seq, siInfo.vcf_ref_seq);

            // STARKA-248 filter invalid indel. TODO: filter this issue earlier (occurs as, e.g. 1D1I which matches ref)
            if (siInfo.vcf_indel_seq == siInfo.vcf_ref_seq) continue;

            static const bool is_use_alt_indel(true);
            _dopt.sicaller_grid().get_somatic_indel(_opt,_dopt,
                                                    normal_sif.sampleOptions,
                                                    tumor_sif.sampleOptions,
                                                    indelKey, indelData, NORMAL,TUMOR,
                                                    is_use_alt_indel,
                                                    siInfo.sindel);

            if (siInfo.sindel.is_output())
            {
                siInfo.indelReportInfo = indelData.getReportInfo();

                // get sample specific info:
                for (unsigned t(0); t<2; ++t)
                {
                    const bool is_include_tier2(t!=0);
                    getAlleleSampleReportInfo(_opt, _dopt, indelKey, normalIndelSampleData, normal_sif.basecallBuffer,
                                              is_include_tier2, is_use_alt_indel,
                                              siInfo.nisri[t]);
                    getAlleleSampleReportInfo(_opt, _dopt, indelKey, tumorIndelSampleData, tumor_sif.basecallBuffer,
                                              is_include_tier2, is_use_alt_indel,
                                              siInfo.tisri[t]);
                }

                pos_t indel_pos(indelKey.pos);
                if (indelKey.type != INDEL::BP_RIGHT)
                {
                    indel_pos -= 1;
                }

                _indelWriter.cacheIndel(indel_pos,siInfo);
                _is_skip_process_pos=false;
            }

#if 0
            /// TODO put this option under runtime control...
            /// TODO setup option so that read keys persist longer when needed for this case...
            ///
            static const bool is_print_indel_evidence(false);

            if (is_print_indel_evidence and is_indel)
            {
                report_os << "INDEL_EVIDENCE " << ik;

                typedef indel_data::score_t::const_iterator siter;
                siter i(id.read_path_lnp.begin()), i_end(id.read_path_lnp.end());
                for (; i!=i_end; ++i)
                {
                    const align_id_t read_id(i->first);
                    const ReadPathScores& lnp(i->second);
                    const ReadPathScores pprob(indel_lnp_to_pprob(_dopt,lnp));
                    const starling_read* srptr(sif.readBuffer.get_read(read_id));

                    report_os << "read key: ";
                    if (nullptr==srptr) report_os << "UNKNOWN_KEY";
                    else            report_os << srptr->key();
                    report_os << "\n"
                              << "read log_lhoods: " << lnp << "\n"
                              << "read pprobs: " << pprob << "\n";
                }
            }
#endif
        }
    }
}



void
strelka_pos_processor::
run_post_call_step(
    const int stage_no,
    const pos_t pos)
{
    if (stage_no != static_cast<int>(_dopt.sfilter.indelRegionStage))
    {
        base_t::run_post_call_step(stage_no, pos);
        return;
    }

    if (! _indelWriter.testPos(pos)) return;

    const LocalRegionStats& was_normal(
        sample(STRELKA_SAMPLE_TYPE::NORMAL).localRegionStatsCollection.getLocalRegionStats(_indelRegionIndexNormal));
    const LocalRegionStats& was_tumor(sample(STRELKA_SAMPLE_TYPE::TUMOR).localRegionStatsCollection.getLocalRegionStats(_indelRegionIndexTumor));

    _indelWriter.addIndelWindowData(_chromName, pos, was_normal, was_tumor, _maxChromDepth);
}


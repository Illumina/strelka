// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Starka
// Copyright (c) 2009-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

///
/// \author Chris Saunders
///

#include "position_somatic_snv.hh"
#include "position_somatic_snv_strand_grid.hh"
#include "position_somatic_snv_strand_grid_vcf.hh"
#include "somatic_indel_grid.hh"
#include "strelka_pos_processor.hh"

#include "blt_util/log.hh"
#include "calibration/scoringmodels.hh"
#include "starling_common/starling_indel_report_info.hh"
#include "starling_common/starling_pos_processor_base_stages.hh"

#include <iomanip>



strelka_pos_processor::
strelka_pos_processor(
    const strelka_options& opt,
    const strelka_deriv_options& dopt,
    const reference_contig_segment& ref,
    const strelka_streams& streams)
    : base_t(opt,dopt,ref,streams,STRELKA_SAMPLE_TYPE::SIZE)
    , _opt(opt)
    , _dopt(dopt)
    , _streams(streams)
    , _scallProcessor(streams.somatic_callable_osptr())
    , _indelWriter(opt, dopt, streams.somatic_indel_osptr())
    , _indelRegionIndexNormal(0)
    , _indelRegionIndexTumor(0)
{
    using namespace STRELKA_SAMPLE_TYPE;

    sample_info& normal_sif(sample(NORMAL));
    sample_info& tumor_sif(sample(TUMOR));

    // set sample-specific parameter overrides:
    normal_sif.sample_opt.min_read_bp_flank = opt.normal_sample_min_read_bp_flank;

    if (opt.is_tumor_sample_min_candidate_indel_reads)
    {
        tumor_sif.sample_opt.min_candidate_indel_reads = opt.tumor_sample_min_candidate_indel_reads;
    }
    if (opt.is_tumor_sample_min_small_candidate_indel_read_frac)
    {
        tumor_sif.sample_opt.min_small_candidate_indel_read_frac = opt.tumor_sample_min_small_candidate_indel_read_frac;
    }

    // setup indel syncronizers:
    {
        double max_candidate_normal_sample_depth(-1.);
        if (dopt.sfilter.is_max_depth())
        {
            if (opt.max_candidate_indel_depth_factor > 0.)
            {
                max_candidate_normal_sample_depth = (opt.max_candidate_indel_depth_factor * dopt.sfilter.max_depth);
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

        indel_sync_data isdata;
        isdata.register_sample(normal_sif.indel_buff,normal_sif.estdepth_buff,normal_sif.estdepth_buff_tier2,
                               normal_sif.sample_opt, max_candidate_normal_sample_depth, NORMAL);
        isdata.register_sample(tumor_sif.indel_buff,tumor_sif.estdepth_buff,tumor_sif.estdepth_buff_tier2,
                               tumor_sif.sample_opt, -1., TUMOR);
        normal_sif.indel_sync_ptr.reset(new indel_synchronizer(opt,ref,isdata,NORMAL));
        tumor_sif.indel_sync_ptr.reset(new indel_synchronizer(opt,ref,isdata,TUMOR));
    }

    // setup indel avg window:
    _indelRegionIndexNormal=normal_sif.wav.add_win(opt.sfilter.indelRegionFlankSize*2);
    _indelRegionIndexTumor=tumor_sif.wav.add_win(opt.sfilter.indelRegionFlankSize*2);
}



void
strelka_pos_processor::
insert_noise_pos(
    const pos_t pos,
    const SiteNoise& sn)
{
    _stageman.validate_new_pos_value(pos,STAGE::READ_BUFFER);
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
    CleanedPileup* normal_cpi_ptr[n_tier] = { &(normal_sif.cpi) , &(_tier2_cpi[NORMAL]) };
    CleanedPileup* tumor_cpi_ptr[n_tier] = { &(tumor_sif.cpi) , &(_tier2_cpi[TUMOR]) };

    for (unsigned t(0); t<n_tier; ++t)
    {
        const bool is_include_tier2(t!=0);
        if (is_include_tier2 && (! _opt.tier2.is_tier2())) continue;
        _pileupCleaner.CleanPileup(normal_sif.bc_buff.get_pos(pos),is_include_tier2,*(normal_cpi_ptr[t]));
        _pileupCleaner.CleanPileup(tumor_sif.bc_buff.get_pos(pos),is_include_tier2,*(tumor_cpi_ptr[t]));
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
        if (_opt.tier2.is_tier2())
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
            _scallProcessor.addToRegion(_chrom_name,output_pos,sgtg);
        }
    }

    // report events:
    //
    bool is_reported_event(false);

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
//        std::ostream& bos(std::cout);

        // have to keep tier1 counts for filtration purposes:
#ifdef SOMATIC_DEBUG
        write_snv_prefix_info_file(_chrom_name,output_pos,ref_base,normald,tumord,log_os);
        log_os << "\n";
#endif

        bos << _chrom_name << '\t'
            << output_pos << '\t'
            << ".";

        static const bool is_write_nqss(false);
        write_vcf_somatic_snv_genotype_strand_grid(_opt,_dopt,sgtg,is_write_nqss,
                                                   *(normal_cpi_ptr[0]),
                                                   *(tumor_cpi_ptr[0]),
                                                   *(normal_cpi_ptr[1]),
                                                   *(tumor_cpi_ptr[1]),
                                                   bos);
        bos << "\n";

        is_reported_event = true;
    }

    if (_opt.is_print_all_site_evidence || (_opt.is_print_evidence && is_reported_event))
    {
        log_os << "TUMOR/NORMAL EVIDENCE pos: " << output_pos << "\n"
               << "is_snv: " << sgtg.is_snv() << "\n"
               << "normal-data:\n" << static_cast<const CleanedPileup&>(*normal_cpi_ptr[0]).cleanedPileup() << "\n"
               << "tumor-data:\n" << static_cast<const CleanedPileup&>(*tumor_cpi_ptr[0]).cleanedPileup() << "\n";
    }
}


void
strelka_pos_processor::
process_pos_variants_impl(const pos_t pos)
{
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

    // with indel synchronization turned on in this model, we should
    // only need to iterate through one sample or the other -- for now
    // we just pick one -- would be nicer in the future to have a way
    // to
    //
    typedef indel_buffer::const_iterator ciiter;
    ciiter i(tumor_sif.indel_buff.pos_iter(pos));
    const ciiter i_end(tumor_sif.indel_buff.pos_iter(pos+1));

    for (; i!=i_end; ++i)
    {
        const indel_key& ik(i->first);

        // don't write breakpoint output:
        if (ik.is_breakpoint()) continue;

        const indel_data& tumor_id(get_indel_data(i));

        if (! tumor_sif.indel_sync().is_candidate_indel(ik,tumor_id)) continue;

        const indel_data* normal_id_ptr(normal_sif.indel_buff.get_indel_data_ptr(ik));
        assert(NULL != normal_id_ptr);
        const indel_data& normal_id(*normal_id_ptr);

        if (normal_id.read_path_lnp.empty() && tumor_id.read_path_lnp.empty()) continue;

        //bool is_indel(false);

        if (_opt.is_somatic_indel())
        {
            // indel_report_info needs to be run first now so that
            // local small repeat info is available to the indel
            // caller

            // get iri from either sample:
            starling_indel_report_info iri;
            get_starling_indel_report_info(ik,tumor_id,_ref,iri);

            // STARKA-248 filter invalid indel. TODO: filter this issue earlier (occurs as, e.g. 1D1I which matches ref)
            if (iri.vcf_indel_seq == iri.vcf_ref_seq) continue;

            double indel_error_prob(0);
            double ref_error_prob(0);
            scoring_models::Instance().get_indel_model().calc_prop(_opt,iri,indel_error_prob,ref_error_prob);

            somatic_indel_call sindel;
            static const bool is_use_alt_indel(true);
            _dopt.sicaller_grid().get_somatic_indel(_opt,_dopt,
                                                    normal_sif.sample_opt,
                                                    tumor_sif.sample_opt,
                                                    indel_error_prob,ref_error_prob,
                                                    ik,normal_id,tumor_id,
                                                    is_use_alt_indel,
                                                    sindel);

            if (sindel.is_output())
            {
                // get sample specific info:
                SomaticIndelVcfInfo siInfo;
                siInfo.sindel = sindel;
                siInfo.iri = iri;

                for (unsigned t(0); t<2; ++t)
                {
                    const bool is_include_tier2(t!=0);
                    get_starling_indel_sample_report_info(_dopt,ik,normal_id,normal_sif.bc_buff,
                                                          is_include_tier2,is_use_alt_indel,
                                                          siInfo.nisri[t]);
                    get_starling_indel_sample_report_info(_dopt,ik,tumor_id,tumor_sif.bc_buff,
                                                          is_include_tier2,is_use_alt_indel,
                                                          siInfo.tisri[t]);
                }

                pos_t indel_pos(ik.pos);
                if (ik.type != INDEL::BP_RIGHT)
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
                    const read_path_scores& lnp(i->second);
                    const read_path_scores pprob(indel_lnp_to_pprob(_dopt,lnp));
                    const starling_read* srptr(sif.read_buff.get_read(read_id));

                    report_os << "read key: ";
                    if (NULL==srptr) report_os << "UNKNOWN_KEY";
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
write_counts(const pos_range& output_report_range) const
{

    std::ostream* report_os_ptr(get_report_osptr());
    if (NULL==report_os_ptr) return;
    std::ostream& report_os(*report_os_ptr);

    for (unsigned i(0); i<STRELKA_SAMPLE_TYPE::SIZE; ++i)
    {
        const sample_info& sif(sample(i));
        const std::string label(STRELKA_SAMPLE_TYPE::get_label(i));

        report_os << std::setprecision(8);
        report_stream_stat(sif.ss,(label+"_ALLSITES_COVERAGE").c_str(),output_report_range,report_os);
        report_stream_stat(sif.used_ss,(label+"_ALLSITES_COVERAGE_USED").c_str(),output_report_range,report_os);

        if (_opt.is_ref_set())
        {
            report_stream_stat(sif.ssn,(label+"_NO_REF_N_COVERAGE").c_str(),output_report_range,report_os);
            report_stream_stat(sif.used_ssn,(label+"_NO_REF_N_COVERAGE_USED").c_str(),output_report_range,report_os);
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

    const win_avg_set& was_normal(sample(STRELKA_SAMPLE_TYPE::NORMAL).wav.get_win_avg_set(_indelRegionIndexNormal));
    const win_avg_set& was_tumor(sample(STRELKA_SAMPLE_TYPE::TUMOR).wav.get_win_avg_set(_indelRegionIndexTumor));

    _indelWriter.addIndelWindowData(pos, was_normal, was_tumor);
}

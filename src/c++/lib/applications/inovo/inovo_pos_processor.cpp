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

#include "extended_pos_data.hh"
#include "inovo_pos_processor.hh"

#include "blt_util/log.hh"
#include "starling_common/starling_indel_error_prob.hh"
#include "starling_common/starling_indel_report_info.hh"
#include "starling_common/starling_pos_processor_base_stages.hh"

#include <iomanip>



inovo_pos_processor::
inovo_pos_processor(
    const inovo_options& opt,
    const inovo_deriv_options& dopt,
    const reference_contig_segment& ref,
    const inovo_streams& streams)
    : base_t(opt,dopt,ref,streams,opt.alignFileOpt.alignmentSampleInfo.size())
    , _opt(opt)
    , _dopt(dopt)
    , _streams(streams)
{
    /// get max proband depth
    double max_candidate_proband_sample_depth(-1.);
    {
        if (dopt.dfilter.is_max_depth())
        {
            if (opt.max_candidate_indel_depth_factor > 0.)
            {
                max_candidate_proband_sample_depth = (opt.max_candidate_indel_depth_factor * dopt.dfilter.max_depth);
            }
        }

        if (opt.max_candidate_indel_depth > 0.)
        {
            if (max_candidate_proband_sample_depth > 0.)
            {
                max_candidate_proband_sample_depth = std::min(max_candidate_proband_sample_depth,static_cast<double>(opt.max_candidate_indel_depth));
            }
            else
            {
                max_candidate_proband_sample_depth = opt.max_candidate_indel_depth;
            }
        }
    }

    using namespace INOVO_SAMPLETYPE;

    // setup indel syncronizers:
    {
        indel_sync_data isdata;
        for (unsigned sampleIndex(0); sampleIndex<_n_samples; ++sampleIndex)
        {
            const bool isProband(_opt.alignFileOpt.alignmentSampleInfo.getSampleInfo(sampleIndex).stype == PROBAND);
            double max_candidate_sample_depth(isProband ? max_candidate_proband_sample_depth : -1);
            sample_info& sif(sample(sampleIndex));
            isdata.register_sample(sif.indel_buff,sif.estdepth_buff,sif.estdepth_buff_tier2,
                                   sif.sample_opt, max_candidate_sample_depth, sampleIndex);
        }

        for (unsigned sampleIndex(0); sampleIndex<_n_samples; ++sampleIndex)
        {
            sample_info& sif(sample(sampleIndex));
            sif.indel_sync_ptr.reset(new indel_synchronizer(opt,ref,isdata,sampleIndex));
        }
    }
}



void
inovo_pos_processor::
process_pos_snp_denovo(const pos_t pos)
{
    using namespace INOVO_SAMPLETYPE;

    const pos_t output_pos(pos+1);
    const char ref_base(_ref.get_base(pos));

    for (unsigned i(0); i<_opt.alignFileOpt.alignmentSampleInfo.size(); ++i)
    {
        const sample_info& sif(sample(i));

    sample_info& normal_sif(sample(NORMAL));
    sample_info& tumor_sif(sample(TUMOR));

    const bool is_dep(_is_dependent_eprob);

    // TODO this is ridiculous -- if the tier2 data scheme works then come back and clean this up:
    static const unsigned n_tier(2);
    std::unique_ptr<extended_pos_data> normald_ptr[n_tier];
    std::unique_ptr<extended_pos_data> tumord_ptr[n_tier];

    extra_position_data* normal_epd_ptr[n_tier] = { &(normal_sif.epd) , &(_tier2_epd[NORMAL]) };
    extra_position_data* tumor_epd_ptr[n_tier] = { &(tumor_sif.epd) , &(_tier2_epd[TUMOR]) };
#if 0
    normald_ptr[0].reset(normal_sif.bc_buff.get_pos(pos),normal_sif.epd,
                         ref_base,_opt,_dpcache,is_dep,is_include_tier2);
    tumord_ptr[0].reset(tumor_sif.bc_buff.get_pos(pos),tumor_sif.epd,
                        ref_base,_opt,_dpcache,is_dep,is_include_tier2);
#endif

    for (unsigned t(0); t<n_tier; ++t)
    {
        const bool is_include_tier2(t!=0);
        if (is_include_tier2 && (! _opt.is_tier2())) continue;
        normald_ptr[t].reset(new extended_pos_data(normal_sif.bc_buff.get_pos(pos),*(normal_epd_ptr[t]),
                                                   ref_base,_opt,_dpcache,is_dep,is_include_tier2));
        tumord_ptr[t].reset(new extended_pos_data(tumor_sif.bc_buff.get_pos(pos),*(tumor_epd_ptr[t]),
                                                  ref_base,_opt,_dpcache,is_dep,is_include_tier2));
    }

#if 0
    sif.ss.update(n_calls);
    sif.used_ss.update(n_used_calls);
    if (pi.ref_base != 'N')
    {
        sif.ssn.update(n_calls);
        sif.used_ssn.update(n_used_calls);
    }
#endif

    denovo_snv_call dsc;

    const extended_pos_info* normal_epi_t2_ptr(NULL);
    const extended_pos_info* tumor_epi_t2_ptr(NULL);
    if (_opt.is_tier2())
    {
        normal_epi_t2_ptr=(&(normald_ptr[1]->good_epi));
        tumor_epi_t2_ptr=(&(tumord_ptr[1]->good_epi));
    }

    position_denovo_snv_call(
        normald_ptr[0]->good_epi,
        tumord_ptr[0]->good_epi,
        normal_epi_t2_ptr,
        tumor_epi_t2_ptr,
        dsc);


    // report events:
    //
    bool is_reported_event(false);

    if (dsc.is_output())
    {
        std::ostream& bos(*_streams.denovo_osptr());
        bos << _chrom_name << '\t'
            << output_pos << '\t'
            << ".";

        write_vcf_denovo_snv(_opt,_dopt,dsc,
                           *(normald_ptr[0]),
                           *(tumord_ptr[0]),
                           *(normald_ptr[1]),
                           *(tumord_ptr[1]),
                           bos);
        bos << "\n";

        is_reported_event = true;
    }
}


void
inovo_pos_processor::
process_pos_variants_impl(const pos_t pos)
{
    try
    {
        process_pos_indel_denovo(pos);
    }
    catch (...)
    {
        log_os << "Exception caught while attempting to call denovo indel at position: " << (pos+1) << "\n";
        throw;
    }

    try
    {
        process_pos_snp_denovo(pos);
    }
    catch (...)
    {
        log_os << "Exception caught while attempting to call denovo SNV at position: " << (pos+1) << "\n";
        throw;
    }

}



void
inovo_pos_processor::
process_pos_indel_denovo(const pos_t pos)
{
#if 0
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

            double indel_error_prob(0);
            double ref_error_prob(0);
            get_indel_error_prob(_opt,iri,indel_error_prob,ref_error_prob);

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
#endif
}



void
inovo_pos_processor::
write_counts(
    const pos_range& output_report_range) const
{
    std::ostream* report_os_ptr(get_report_osptr());
    if (nullptr==report_os_ptr) return;
    std::ostream& report_os(*report_os_ptr);

    for (unsigned i(0); i<_opt.alignFileOpt.alignmentSampleInfo.size(); ++i)
    {
        const sample_info& sif(sample(i));
        const std::string label(INOVO_SAMPLETYPE::get_label(i));

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

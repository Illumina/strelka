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

#include "inovo_pos_processor.hh"

#include "blt_util/log.hh"
#include "starling_common/starling_indel_error_prob.hh"
#include "starling_common/starling_indel_report_info.hh"

#include <iomanip>
#include "denovo_snv_caller.hh"
#include "denovo_snv_call_vcf.hh"



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

    // skip site if proband depth is zero:
    {
        const unsigned probandIndex(_opt.alignFileOpt.alignmentSampleInfo.getTypeIndexList(PROBAND)[0]);
        const CleanedPileup& probandCpi(sample(probandIndex).cpi);

        // note this is a more expansive skipping criteria then we use for germline calling
        // (this is because there's no gvcf output)
        if (probandCpi.cleanedPileup().calls.empty()) return;
    }

    // finish formatting pileups for all samples:
    for (unsigned i(0); i<_n_samples; ++i)
    {
        sample_info& sif(sample(i));
        _pileupCleaner.CleanPileupErrorProb(sif.cpi);
    }

    const pos_t output_pos(pos+1);
    //const char ref_base(_ref.get_base(pos));

    // make cleaned pileups of all samples easily accessible to the variant caller:
    std::vector<const CleanedPileup*> pileups(_n_samples);
    for (unsigned i(0); i<_n_samples; ++i)
    {
        pileups[i] = &(sample[i]);
    }

    denovo_snv_call dsc;

    get_denovo_snv_call(
        _opt,
        _opt.alignFileOpt.alignmentSampleInfo,
        pileups,
        dsc);

    // report events:
    //
    if (dsc.is_output())
    {
        std::ostream& bos(*_streams.denovo_osptr());
        bos << _chrom_name << '\t'
            << output_pos << '\t'
            << ".";

        denovo_snv_call_vcf(
            _opt,_dopt,
            _opt.alignFileOpt.alignmentSampleInfo,
            pileups,
            dsc,
            bos);
        bos << "\n";
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

    typedef indel_buffer::const_iterator ciiter;
    ciiter i(proband_sif.indel_buff.pos_iter(pos));
    const ciiter i_end(proband_sif.indel_buff.pos_iter(pos+1));

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

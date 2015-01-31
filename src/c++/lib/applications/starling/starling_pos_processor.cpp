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

#include "starling_pos_processor.hh"
#include "blt_common/position_nonref_test.hh"
#include "blt_common/position_nonref_2allele_test.hh"
//#include "blt_common/position_snp_call_lrt.hh"
#include "blt_common/ref_context.hh"

#include <iomanip>



static
void
report_counts(const snp_pos_info& pi,
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



static
void
write_bsnp_diploid_allele(const blt_options& client_opt,
                          const blt_streams& client_io,
                          const std::string& seq_name,
                          const pos_t output_pos,
                          const char ref,
                          const unsigned n_used_calls,
                          const unsigned n_unused_calls,
                          const snp_pos_info& good_pi,
                          const diploid_genotype& dgt,
                          const unsigned hpol = 0)
{
    std::ostream& os(*client_io.bsnp_diploid_allele_osptr());

    write_snp_prefix_info_file(seq_name,output_pos,ref,n_used_calls,n_unused_calls,os);
    os << "\t";
    write_diploid_genotype_allele(client_opt,good_pi,dgt,os,hpol);
    os << "\n";
}


starling_pos_processor::
starling_pos_processor(
    const starling_options& opt,
    const starling_deriv_options& dopt,
    const reference_contig_segment& ref,
    const starling_streams_base& client_io)
    : base_t(opt,dopt,ref,client_io,1)
{
    // setup indel syncronizers:
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

        indel_sync_data isdata;
        isdata.register_sample(normal_sif.indel_buff,normal_sif.estdepth_buff,normal_sif.estdepth_buff_tier2,
                               normal_sif.sample_opt, max_candidate_normal_sample_depth, 0);
        normal_sif.indel_sync_ptr.reset(new indel_synchronizer(opt, ref, isdata, 0));
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
        process_pos_snp_single_sample_impl(pos,sample_no);
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

    const snp_pos_info& pi(sif.bc_buff.get_pos(pos));
    const snp_pos_info& good_pi(sif.epd.good_pi);
    const pos_t output_pos(pos+1);

    // note multi-sample status -- can still be called only for one sample
    // and only for sample 0. working on generalization:
    //
    if (sample_no!=0) return;

    if (pi.calls.empty()) return;

    adjust_joint_eprob(_client_opt,_dpcache,good_pi,sif.epd.dependent_eprob,_is_dependent_eprob);

    const extended_pos_info good_epi(good_pi,sif.epd.dependent_eprob);

    // get fraction of filtered bases:
#if 0
    const double filter_fraction(static_cast<double>(n_unused_calls)/static_cast<double>(n_calls));
    const bool is_overfilter(filter_fraction > _client_opt.max_basecall_filter_fraction);
#endif

    // delay writing any snpcalls so that anomaly tests can (optionally) be applied as filters:
    //
    nonref_test_call nrc;
    //lrt_snp_call lsc;
    _site_info.dgt.reset();
    //monoploid_genotype mgt;
    //std::unique_ptr<nploid_genotype> ngt_ptr;

    // check whether we're in a haploid region:
    _site_info.dgt.ploidy=(get_ploidy(pos));

    if (_client_opt.is_counts)
    {
        report_counts(good_pi,_site_info.n_unused_calls,output_pos,*_client_io.counts_osptr());
    }

    if (_client_opt.is_nonref_test() || _client_opt.is_nonref_sites())
    {
        position_nonref_2allele_test(good_pi,
                                     _client_opt,
                                     _client_opt.is_nonref_sites(),
                                     nrc);
#if 0
        static const bool is_mle_freq(false);

        position_nonref_test(good_pi,
                             _client_opt.nonref_variant_rate,
                             _client_opt.min_nonref_freq,
                             is_mle_freq,
                             nrc);
#endif

    }

#if 0
    if (_client_opt.is_lsnp)
    {
        position_snp_call_lrt(_client_opt.lsnp_alpha,good_pi,lsc);
    }
#endif
    if (_client_opt.is_bsnp_diploid())
    {
        _client_dopt.pdcaller().position_snp_call_pprob_digt(
            _client_opt,good_epi,_site_info.dgt,_client_opt.is_all_sites());
    }
#if 0
    if (_client_opt.is_bsnp_monoploid)
    {
        position_snp_call_pprob_monogt(_client_opt.bsnp_monoploid_theta,good_pi,mgt);
    }
    if (_client_opt.is_bsnp_nploid)
    {
        ngt_ptr.reset(new nploid_genotype(*_ninfo));
        position_snp_call_pprob_nploid(_client_opt.bsnp_nploid_snp_prob,good_pi,*_ninfo,*ngt_ptr);
    }
#endif

    //    const bool is_snp(nrc.is_snp || lsc.is_snp || _site_info.dgt.is_snp || mgt.is_snp || (ngt_ptr.get() && ngt_ptr->is_snp));
    const bool is_snp(nrc.is_snp || _site_info.dgt.is_snp);

    // find anomalies:
    //
#if 0
    bool is_pos_adis(false);
    bool is_pos_acov(false);

    if ((_client_opt.is_adis_table || _client_opt.is_adis_lrt) && is_snp)
    {
        if (_client_opt.is_adis_table)
        {
            is_pos_adis = (is_pos_adis || position_strand_distro_anomaly(_client_opt.adis_table_alpha,good_pi,_ws));
        }
        if (_client_opt.is_adis_lrt)
        {
            is_pos_adis = (is_pos_adis || position_strand_distro_anomaly_lrt(_client_opt.adis_lrt_alpha,good_pi));
        }
    }
    if (_client_opt.is_acov)
    {
        is_pos_acov = position_strand_coverage_anomaly(_client_opt.acov_alpha,pi);
    }
#endif

    //const bool is_anomaly(is_pos_adis || is_pos_acov);
    //const bool is_filter_snp(is_overfilter || (_client_opt.is_filter_anom_calls && is_anomaly));

    //    const bool is_nf_snp(is_snp && (! is_filter_snp));
    if (is_snp)
    {
        if (_client_opt.is_compute_hapscore)
        {
            _site_info.hapscore=get_hapscore(pi.hap_set);
        }

        // do calculate VQSR metrics
        if (_client_opt.is_compute_germline_VQSRmetrics())
        {
            _site_info.MQ               = pi.get_rms_mq();
            _site_info.ReadPosRankSum   = pi.get_read_pos_ranksum();
            _site_info.MQRankSum        = pi.get_mq_ranksum();
            _site_info.BaseQRankSum     = pi.get_baseq_ranksum();
            _site_info.rawPos           = pi.get_raw_pos();
            _site_info.avgBaseQ         = pi.get_raw_baseQ();
        }

        // hpol filter
        _site_info.hpol=get_snp_hpol_size(pos,_ref);
    }

    if (_client_opt.is_all_sites())
    {
#if 0
        const diploid_genotype* dgt_ptr(&_site_info.dgt);
        if (is_filter_snp)
        {
            dgt_ptr=&get_empty_dgt(pi.ref_base);
        }
#endif

        //Add site to gvcf
        if (_client_opt.gvcf.is_gvcf_output())
        {
            _site_info.init(pos,pi.get_ref_base(),good_pi,_client_opt.used_allele_count_min_qscore);
            _gvcfer->add_site(_site_info);
        }


        if (_client_opt.is_bsnp_diploid_allele_file)
        {
            write_bsnp_diploid_allele(_client_opt,_client_io,_chrom_name,output_pos,pi.get_ref_base(),_site_info.n_used_calls,_site_info.n_unused_calls,good_pi,_site_info.dgt,_site_info.hpol);
        }
    }

    if (_client_opt.is_nonref_sites())
    {
        std::ostream& bos(*_client_io.nonref_sites_osptr());
        write_snp_prefix_info_file(_chrom_name,output_pos,pi.get_ref_base(),_site_info.n_used_calls,_site_info.n_unused_calls,bos);
        bos << "\t";
        write_nonref_2allele_test(_client_opt,good_pi,nrc,bos);
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
            std::ostream& bos(*_client_io.nonref_test_osptr());
            write_snp_prefix_info_file(_chrom_name,output_pos,pi.get_ref_base(),_site_info.n_used_calls,_site_info.n_unused_calls,bos);
            bos << "\t";
            write_nonref_2allele_test(_client_opt,good_pi,nrc,bos);
#if 0
            write_nonref_test(_client_opt,good_pi,nrc,bos);
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
        if (_site_info.dgt.is_snp)
        {
            if (_client_opt.is_bsnp_diploid_file)
            {
                std::ostream& bos(*_client_io.bsnp_diploid_osptr());
                write_snp_prefix_info_file(_chrom_name,output_pos,pi.get_ref_base(),_site_info.n_used_calls,_site_info.n_unused_calls,bos);
                bos << "\t";
                write_diploid_genotype_snp(_client_opt,good_pi,_site_info.dgt,bos,_site_info.hpol);
                bos << "\n";
            }

            // this needs to be updated no matter where the snp-call is written to:
            if (_is_variant_windows) _variant_print_pos.insert(pos);
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

#if 0
    if (is_anomaly && (! _client_opt.is_filter_anom_calls))
    {
        if (is_pos_adis) report_os << "ANOM_DIS pos: " << output_pos << "\n";
        if (is_pos_acov) report_os << "ANOM_COV pos: " << output_pos << "\n";

        is_reported_event = true;
    }
#endif

    if (_client_opt.is_print_all_site_evidence || (_client_opt.is_print_evidence && is_reported_event))
    {
        report_os << "EVIDENCE pos: " << output_pos << "\n"
                  << "is_snp: " << is_snp << "\n"
                  << pi << "\n";
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

    if (_client_opt.is_ref_set())
    {
        report_stream_stat(sif.ssn,"NO_REF_N_COVERAGE",output_report_range,report_os);
        report_stream_stat(sif.used_ssn,"NO_REF_N_COVERAGE_USED",output_report_range,report_os);
    }
}

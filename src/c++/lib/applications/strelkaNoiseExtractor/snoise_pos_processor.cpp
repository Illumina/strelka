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

#include "snoise_pos_processor.hh"

#include <algorithm>
#include <array>
#include <iomanip>



void
snoise_pos_processor::
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



void
snoise_pos_processor::
process_pos_snp_snoise(
    const pos_t pos)
{
    const unsigned sample_no(0);
    sample_info& sif(sample(sample_no));

    snp_pos_info null_pi;
    snp_pos_info* pi_ptr(sif.bc_buff.get_pos(pos));
    if (NULL==pi_ptr) pi_ptr=&null_pi;
    snp_pos_info& pi(*pi_ptr);

    const unsigned n_calls(pi.calls.size());
    const unsigned n_spandel(pi.n_spandel);
    const unsigned n_submapped(pi.n_submapped);

    const pos_t output_pos(pos+1);

    pi.set_ref_base(_ref.get_base(pos));

    std::ostream& bos(*_client_io.snoise_osptr());

    // for all but coverage-tests, we use a high-quality subset of the basecalls:
    //
    snp_pos_info& good_pi(sif.epd.good_pi);
    good_pi.clear();
    good_pi.set_ref_base(pi.get_ref_base());
    for (unsigned i(0); i<n_calls; ++i)
    {
        if (pi.calls[i].is_call_filter) continue;
        good_pi.calls.push_back(pi.calls[i]);
    }

    const unsigned n_used_calls=(good_pi.calls.size());
    const unsigned n_unused_calls=(n_calls-_site_info.n_used_calls);


    // note multi-sample status -- can still be called only for one sample
    // and only for sample 0. working on generalization:
    //
    if (sample_no!=0) return;

    if (pi.calls.empty()) return;

    // approximate the the site is not germline and contains noise instead:
    if (n_used_calls < 15) return;

    std::array<unsigned,N_BASE> base_count;
    std::fill(base_count.begin(),base_count.end(),0);

    for (const auto& bc : good_pi.calls)
    {
        assert(bc.base_id!=BASE_ID::ANY);
        base_count[bc.base_id]++;
    }


    report_counts(good_pi,_site_info.n_unused_calls,output_pos,*_client_io.counts_osptr());


    static
    void
    report_counts(const snp_pos_info& pi,
                  const unsigned n_unused_calls,
                  const pos_t output_pos,
                  std::ostream& os)
    {

        for (unsigned i(0); i<N_BASE; ++i) base_count[i] = 0;



        os << output_pos << '\t';
        for (unsigned i(0); i<N_BASE; ++i)
        {
            os << base_count[i] << '\t';
        }
        os << n_unused_calls << '\n';
    }


#if 0
    // for all but coverage-tests, we use a high-quality subset of the basecalls:
    //
    snp_pos_info& good_pi(sif.epd.good_pi);
    good_pi.clear();
    good_pi.set_ref_base(pi.get_ref_base());
    for (unsigned i(0); i<n_calls; ++i)
    {
        if (pi.calls[i].is_call_filter) continue;
        good_pi.calls.push_back(pi.calls[i]);
    }

    _site_info.n_used_calls=(good_pi.calls.size());
    _site_info.n_unused_calls=(n_calls-_site_info.n_used_calls);

    sif.ss.update(n_calls);
    sif.used_ss.update(_site_info.n_used_calls);
    if (pi.get_ref_base() != 'N')
    {
        sif.ssn.update(n_calls);
        sif.used_ssn.update(_site_info.n_used_calls);
        sif.wav.insert(pos,_site_info.n_used_calls,_site_info.n_unused_calls,n_spandel,n_submapped);
    }
    else
    {
        sif.wav.insert_null(pos);
    }

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
        _client_dopt.pdcaller().position_snp_call_pprob_digt(_client_opt,good_epi,_site_info.dgt,_client_opt.is_all_sites());
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
        if (_client_opt.is_compute_VQSRmetrics || _client_opt.calibration_model!="default")
        {
            _site_info.MQ               = pi.get_rms_mq();
            _site_info.ReadPosRankSum   = pi.get_read_pos_ranksum();
            _site_info.MQRankSum        = pi.get_mq_ranksum();
            _site_info.BaseQRankSum     = pi.get_baseq_ranksum();
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
        if (_client_opt.is_gvcf_output())
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
#endif
}

// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/downloads/sequencing/licenses/>.
//

/// \file
///
/// \author Chris Saunders
///

#include "extended_pos_data.hh"
#include "position_somatic_snv.hh"
#include "position_somatic_snv_grid.hh"
#include "position_somatic_snv_strand_grid.hh"
#include "somatic_indel.hh"
#include "somatic_indel_grid.hh"
#include "strelka_pos_processor.hh"

#include "blt_util/log.hh"
#include "starling_common/starling_indel_error_prob.hh"
#include "starling_common/starling_indel_report_info.hh"

#include <iomanip>



strelka_pos_processor::
strelka_pos_processor(const strelka_options& opt,
                      const strelka_deriv_options& dopt,
                      const reference_contig_segment& ref,
                      const strelka_streams& client_io)
    : base_t(opt,dopt,ref,client_io,STRELKA_SAMPLE_TYPE::SIZE)
    , _opt(opt)
    , _dopt(dopt)
    , _client_io(client_io) {

    using namespace STRELKA_SAMPLE_TYPE;

    sample_info& normal_sif(sample(NORMAL));
    sample_info& tumor_sif(sample(TUMOR));

    // set sample-specific parameter overrides:
    normal_sif.sample_opt.min_read_bp_flank = opt.normal_sample_min_read_bp_flank;

    if (opt.is_tumor_sample_min_candidate_indel_reads) {
        tumor_sif.sample_opt.min_candidate_indel_reads = opt.tumor_sample_min_candidate_indel_reads;
    }
    if (opt.is_tumor_sample_min_small_candidate_indel_read_frac) {
        tumor_sif.sample_opt.min_small_candidate_indel_read_frac = opt.tumor_sample_min_small_candidate_indel_read_frac;
    }

    // setup indel syncronizers:
    indel_sync_data isdata;
    isdata.register_sample(normal_sif.indel_buff,normal_sif.estdepth_buff,normal_sif.sample_opt,NORMAL);
    isdata.register_sample(tumor_sif.indel_buff,tumor_sif.estdepth_buff,tumor_sif.sample_opt,TUMOR);
    normal_sif.indel_sync_ptr = (new indel_synchronizer(isdata,NORMAL));
    tumor_sif.indel_sync_ptr = (new indel_synchronizer(isdata,TUMOR));
}



strelka_pos_processor::
~strelka_pos_processor() {
    for (unsigned i(0); i<_n_samples; ++i) {
        delete sample(i).indel_sync_ptr;
    }
}




void
strelka_pos_processor::
process_pos_snp_somatic(const pos_t pos) {

    using namespace STRELKA_SAMPLE_TYPE;

    const pos_t output_pos(pos+1);
    const char ref_base(_ref.get_base(pos));

    sample_info& normal_sif(sample(NORMAL));
    sample_info& tumor_sif(sample(TUMOR));

    const bool is_dep(_is_dependent_eprob);

    // TODO this is ridiculous -- if the tier2 data scheme works then come back and clean this up:
    static const unsigned n_tier(2);
    std::auto_ptr<extended_pos_data> normald_ptr[n_tier];
    std::auto_ptr<extended_pos_data> tumord_ptr[n_tier];

    extra_position_data* normal_epd_ptr[n_tier] = { &(normal_sif.epd) , &(_tier2_epd[NORMAL]) };
    extra_position_data* tumor_epd_ptr[n_tier] = { &(tumor_sif.epd) , &(_tier2_epd[TUMOR]) };
#if 0
    normald_ptr[0].reset(normal_sif.bc_buff.get_pos(pos),normal_sif.epd,
                         ref_base,_opt,_dpcache,is_dep,is_include_tier2);
    tumord_ptr[0].reset(tumor_sif.bc_buff.get_pos(pos),tumor_sif.epd,
                        ref_base,_opt,_dpcache,is_dep,is_include_tier2);
#endif

    for (unsigned t(0); t<n_tier; ++t) {
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
    if (pi.ref_base != 'N') {
        sif.ssn.update(n_calls);
        sif.used_ssn.update(n_used_calls);
    }
#endif

    // note single-sample anomaly filtration won't apply here (more of
    // a vestigal blt feature anyway)
    //

    // retain original blt loop structure from the single-sample case
    // to allow for multiple interacting tests at one site
    //

    //    somatic_snv_genotype sgt;
    somatic_snv_genotype_grid sgtg;

#if 0
    if (_opt.is_somatic_snv()) {
        const extended_pos_info* normal_epi_t2_ptr(NULL);
        const extended_pos_info* tumor_epi_t2_ptr(NULL);
        if (_opt.is_tier2()) {
            normal_epi_t2_ptr=(&(normald_t2_ptr->good_epi));
            tumor_epi_t2_ptr=(&(tumord_t2_ptr->good_epi));
        }
        _dopt.sscaller().position_somatic_snv_call(normald.good_epi,
                                                   tumord.good_epi,
                                                   normal_epi_t2_ptr,
                                                   tumor_epi_t2_ptr,
                                                   sgt);
    }
#else
#if 0
    if (_opt.is_somatic_snv()) {
        const extended_pos_info* normal_epi_t2_ptr(NULL);
        const extended_pos_info* tumor_epi_t2_ptr(NULL);
        if (_opt.is_tier2()) {
            normal_epi_t2_ptr=(&(normald_t2_ptr->good_epi));
            tumor_epi_t2_ptr=(&(tumord_t2_ptr->good_epi));
        }
        _dopt.sscaller_grid().position_somatic_snv_call(normald.good_epi,
                                                        tumord.good_epi,
                                                        normal_epi_t2_ptr,
                                                        tumor_epi_t2_ptr,
                                                        sgtg);
    }
#else
    if (_opt.is_somatic_snv()) {
        const extended_pos_info* normal_epi_t2_ptr(NULL);
        const extended_pos_info* tumor_epi_t2_ptr(NULL);
        if (_opt.is_tier2()) {
            normal_epi_t2_ptr=(&(normald_ptr[1]->good_epi));
            tumor_epi_t2_ptr=(&(tumord_ptr[1]->good_epi));
        }
        _dopt.sscaller_strand_grid().position_somatic_snv_call(normald_ptr[0]->good_epi,
                                                               tumord_ptr[0]->good_epi,
                                                               normal_epi_t2_ptr,
                                                               tumor_epi_t2_ptr,
                                                               sgtg);
    }
#endif
#endif

    const bool is_snv(sgtg.is_snv);

    // report events:
    //
    bool is_reported_event(false);

    if (is_snv) {
#if 0
        if (sgt.is_snv) {
            std::ostream& bos(*_client_io.somatic_snv_osptr());
            const extended_pos_data& ndata( (sgt.tier==0) ? normald : *normald_t2_ptr );
            const extended_pos_data& tdata( (sgt.tier==0) ? tumord : *tumord_t2_ptr );
            //            write_snv_prefix_info_file(_chrom_name,output_pos,ref_base,ndata,tdata,bos);

            // have to keep tier1 counts for filtration purposes:
#ifdef SOMATIC_DEBUG
            write_snv_prefix_info_file(_chrom_name,output_pos,ref_base,normald,tumord,log_os);
            log_os << "\n";
#endif
            write_snv_prefix_info_file(_chrom_name,output_pos,ref_base,normald,tumord,bos);
            bos << "\t";
            write_somatic_snv_genotype(_opt,sgt,ndata.epd.good_pi,tdata.epd.good_pi,bos);
            bos << "\n";
        }
#endif
        if (sgtg.is_snv) {
            std::ostream& bos(*_client_io.somatic_snv_osptr());
            //            const extended_pos_data& ndata( (sgtg.tier==0) ? normald : *normald_t2_ptr );
            //            const extended_pos_data& tdata( (sgtg.tier==0) ? tumord : *tumord_t2_ptr );
            //            write_snv_prefix_info_file(_chrom_name,output_pos,ref_base,ndata,tdata,bos);

            // have to keep tier1 counts for filtration purposes:
#ifdef SOMATIC_DEBUG
            write_snv_prefix_info_file(_chrom_name,output_pos,ref_base,normald,tumord,log_os);
            log_os << "\n";
#endif

            bos << _chrom_name << '\t'
                << output_pos << '\t'
                << ".";

            write_vcf_somatic_snv_genotype_strand_grid(_opt,sgtg,
                                                       *(normald_ptr[0]),
                                                       *(tumord_ptr[0]),
                                                       *(normald_ptr[1]),
                                                       *(tumord_ptr[1]),
                                                       bos);
            bos << "\n";
        }
        is_reported_event = true;
    }

    if (_opt.is_print_all_site_evidence || (_opt.is_print_evidence && is_reported_event)) {
        log_os << "TUMOR/NORMAL EVIDENCE pos: " << output_pos << "\n"
               << "is_snv: " << is_snv << "\n"
               << "normal-data:\n" << normald_ptr[0]->epd.good_pi << "\n"
               << "tumor-data:\n" << tumord_ptr[0]->epd.good_pi << "\n";
    }
}



void
strelka_pos_processor::
process_pos_indel_somatic(const pos_t pos) {

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

    for (; i!=i_end; ++i) {
        const indel_key& ik(i->first);

        // don't write breakpoint output:
        if(ik.is_breakpoint()) continue;

        const indel_data& tumor_id(get_indel_data(i));
        if (! tumor_sif.indel_sync().is_candidate_indel(_opt,ik,tumor_id)) continue;

        const indel_data* normal_id_ptr(normal_sif.indel_buff.get_indel_data_ptr(ik));
        assert(NULL != normal_id_ptr);
        const indel_data& normal_id(*normal_id_ptr);

        if (normal_id.read_path_lnp.empty() && tumor_id.read_path_lnp.empty()) continue;

        //bool is_indel(false);

        if (_opt.is_somatic_indel()) {
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
#ifdef USE_ORIG_SINDEL
            static const bool is_use_alt_indel(true);
            _dopt.sicaller().get_somatic_indel(_opt,_dopt,
                                               indel_error_prob,ref_error_prob,
                                               ik,normal_id,tumor_id,
                                               is_use_alt_indel,
                                               sindel);
#else
            static const bool is_use_alt_indel(true);
            _dopt.sicaller_grid().get_somatic_indel(_opt,_dopt,
                                                    normal_sif.sample_opt,
                                                    tumor_sif.sample_opt,
                                                    indel_error_prob,ref_error_prob,
                                                    ik,normal_id,tumor_id,
                                                    is_use_alt_indel,
                                                    sindel);
#endif

            if (sindel.is_indel) {
                //is_indel=true;

                // get sample specific info:
                starling_indel_sample_report_info normal_isri[2];
                starling_indel_sample_report_info tumor_isri[2];
                for (unsigned t(0); t<2; ++t) {
                    const bool is_include_tier2(t!=0);
                    get_starling_indel_sample_report_info(_dopt,ik,normal_id,normal_sif.bc_buff,
                                                          is_include_tier2,is_use_alt_indel,
                                                          normal_isri[t]);
                    get_starling_indel_sample_report_info(_dopt,ik,tumor_id,tumor_sif.bc_buff,
                                                          is_include_tier2,is_use_alt_indel,
                                                          tumor_isri[t]);
                }

                pos_t indel_pos(ik.pos);
                if (ik.type != INDEL::BP_RIGHT) {
                    indel_pos -= 1;
                }

                const pos_t output_pos(indel_pos+1);

                std::ostream& bos(*_client_io.somatic_indel_osptr());
                bos << _chrom_name << '\t'
                    << output_pos << '\t' << '.';
                write_somatic_indel_vcf_grid(sindel,iri,normal_isri,tumor_isri,bos);
                bos << '\n';

                _variant_print_pos.insert(indel_pos);
            }

#if 0
            /// TODO put this option under runtime control...
            /// TODO setup option so that read keys persist longer when needed for this case...
            ///
            static const bool is_print_indel_evidence(false);

            if (is_print_indel_evidence and is_indel) {
                report_os << "INDEL_EVIDENCE " << ik;

                typedef indel_data::score_t::const_iterator siter;
                siter i(id.read_path_lnp.begin()), i_end(id.read_path_lnp.end());
                for (; i!=i_end; ++i) {
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
write_counts(const pos_range& output_report_range) const {

    std::ostream* report_os_ptr(get_report_osptr());
    if (NULL==report_os_ptr) return;
    std::ostream& report_os(*report_os_ptr);

    for (unsigned i(0); i<STRELKA_SAMPLE_TYPE::SIZE; ++i) {
        const sample_info& sif(sample(i));
        const std::string label(STRELKA_SAMPLE_TYPE::get_label(i));

        report_os << std::setprecision(8);
        report_stream_stat(sif.ss,(label+"_ALLSITES_COVERAGE").c_str(),output_report_range,report_os);
        report_stream_stat(sif.used_ss,(label+"_ALLSITES_COVERAGE_USED").c_str(),output_report_range,report_os);

        if (_opt.is_ref_set()) {
            report_stream_stat(sif.ssn,(label+"_NO_REF_N_COVERAGE").c_str(),output_report_range,report_os);
            report_stream_stat(sif.used_ssn,(label+"_NO_REF_N_COVERAGE_USED").c_str(),output_report_range,report_os);
        }
    }
}

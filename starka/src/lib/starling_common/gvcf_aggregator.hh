// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2012 Illumina, Inc.
//
// This software is covered by the "Illumina Genome Analyzer Software
// License Agreement" and the "Illumina Source Code License Agreement",
// and certain third party copyright/licenses, and any user of this
// source file is bound by the terms therein (see accompanying files
// Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
// Illumina_Source_Code_License_Agreement.pdf and third party
// copyright/license notices).
//
//

/// \file
///
/// \author Chris Saunders
///

#ifndef __GVCF_AGGREGATOR_HH
#define __GVCF_AGGREGATOR_HH


#include "blt_common/position_snp_call_pprob_digt.hh"
#include "starling_common/starling_indel_call_pprob_digt.hh"
#include "starling_common/starling_shared.hh"

#include <iostream>



struct indel_info {

    void
    init(const pos_t init_pos,
         const indel_key& init_ik,
         const starling_diploid_indel_core& init_dindel,
         const starling_indel_report_info& init_iri,
         const starling_indel_sample_report_info& init_isri)
    { 
        pos=(init_pos);
        ik=(init_ik);
        dindel=(init_dindel);
        iri=(init_iri);
        isri=(init_isri);
    }

    pos_t pos;
    indel_key ik;
    starling_diploid_indel_core dindel;
    starling_indel_report_info iri;
    starling_indel_sample_report_info isri;
};



struct site_info {
    void
    init(const pos_t init_pos,
         const char init_ref,
         const unsigned init_n_used_calls,
         const unsigned init_n_unused_calls,
         const snp_pos_info& init_good_pi,
         const diploid_genotype& init_dgt,
         const bool init_is_nf_snp,
         const double init_sb,
         const unsigned init_hpo) {

        pos=(init_pos);
        ref=(init_ref);
        n_used_calls=(init_n_used_calls);
        n_unused_calls=(init_n_unused_calls);
        good_pi=init_good_pi;
        dgt=init_dgt;
        is_nf_snp=init_is_nf_snp;
        sb=init_sb;
        hpo=init_hpo;
    }

    pos_t pos;
    char ref;
    unsigned n_used_calls;
    unsigned n_unused_calls;
    snp_pos_info good_pi;
    diploid_genotype dgt;
    bool is_nf_snp;
    double sb;
    unsigned hpo;
};



///
/// Assembles all site and indel call information into a consistent set, blocks output
/// and writes to a VCF stream
///

struct gvcf_aggregator {

    gvcf_aggregator(const starling_options& opt,
                    const pos_range& report_range,
                    std::ostream* os);

    ~gvcf_aggregator();

    void
    set_chrom_name(const char* chrom) {
        _chrom=chrom;
    }

    void
    add_site(const pos_t pos,
             const char ref,
             const unsigned n_used_calls,
             const unsigned n_unused_calls,
             const snp_pos_info& good_pi,
             const diploid_genotype& dgt,
             const bool is_nf_snp,
             const double sb,
             const unsigned hpol);

    void
    add_indel(const pos_t pos,
              const indel_key ik,
              const starling_diploid_indel_core& dindel,
              const starling_indel_report_info& iri,
              const starling_indel_sample_report_info& isri);

private:

    // resolve a set of overlapping indel and site calls:
    void process_overlaps();
    

    // initial policy is to write nothing at empty sites. why?
    //
    // (1) gatk does it
    // (2) if gVCF output is every turned off, the output would be ridiculous -- maybe this should be a gVCf only thing?
    //
    //void
    //write_empty_pos();

    //void
    //update_pos();

    const starling_options& _opt;
    const char* _chrom;

    const known_pos_range _report_range;

    //pos_t _head_pos;  // all input has been recieved up to this pos
    //pos_t _write_pos; // output has been written out to this pos

    pos_t _indel_end_pos;

    unsigned _indel_buffer_size;
    std::vector<indel_info> _indel_buffer;

    unsigned _site_buffer_size;
    std::vector<site_info> _site_buffer;

    // convenient reference to gvcf stream from opt:
    std::ostream* _osptr;
};


#endif

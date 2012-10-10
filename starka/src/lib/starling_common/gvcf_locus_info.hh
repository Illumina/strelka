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

#ifndef __GVCF_LOCUS_INFO_HH
#define __GVCF_LOCUS_INFO_HH


#include "blt_common/position_snp_call_pprob_digt.hh"
#include "starling_common/starling_indel_call_pprob_digt.hh"
#include "starling_common/starling_shared.hh"

#include <bitset>


namespace VCF_FILTERS {

    enum index_t {
        HighDepth,
        LowGQX,
        HighSB,
        HighHpol,
        HighBaseFilt,
        HighRepeatCount,
        IndelConflict,
        SiteConflict,
        SIZE
    };

    inline 
    const char*
    get_label(const index_t i) {
        switch(i) {
        case HighDepth: return "HighDepth";
        case LowGQX: return "LowGQX";
        case HighSB: return "HighSB";
        case HighHpol: return "HighHpol";
        case HighBaseFilt: return "HighBaseFilt";
        case HighRepeatCount: return "HighRepeatCount";
        case IndelConflict: return "IndelConflict";
        case SiteConflict: return "SiteConflict";
        default:
            assert(0);
            return NULL;
        }
    }
}



struct shared_modifiers {
    
    shared_modifiers()
    {}

    void
    set_filter(const VCF_FILTERS::index_t i) {
        filters[i] = 1;
    }

    std::bitset<VCF_FILTERS::SIZE> filters;
};


struct indel_modifiers : public shared_modifiers {

};

struct site_modifiers : public shared_modifiers {

};



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

    indel_modifiers imod;
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

    site_modifiers smod;
};



#endif

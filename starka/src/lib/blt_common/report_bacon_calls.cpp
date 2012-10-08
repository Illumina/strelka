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

#include "blt_common/report_bacon_calls.hh"

#include "blt_util/log.hh"

#include <cassert>
#include <cmath>
#include <cstdlib>

#include <iomanip>
#include <iostream>



const double bacon_max2_score_print_thresh(2);



static
bacon_snp_t
bacon_snp_class(const unsigned max_base,
                const unsigned max2_base,
                const unsigned ref_base,
                const bool is_het){

    if(max_base == ref_base){
        return BACON_SNP_HET1;
    } else {
        if(is_het){
            if(max2_base == ref_base){
                return BACON_SNP_HET2;
            } else {
                return BACON_SNP_HETO;
            }
        } else {
            return BACON_SNP_DIFF;
        }
    }
}



static
const char*
bacon_snp_label(const bacon_snp_t s){

    switch(s) {
    case BACON_SNP_DIFF: return "SNP_diff";
    case BACON_SNP_HET1: return "SNP_het1";
    case BACON_SNP_HET2: return "SNP_het2";
    case BACON_SNP_HETO: return "SNP_other_het";
    default:
        log_os << "ERROR:: bacon_snp_label - unknown bacon snp type\n";
        exit(EXIT_FAILURE);
    }
}



void
get_bacon_scores(const blt_options& client_opt,
                 const snp_pos_info& pi,
                 const unsigned unused_count,
                 bacon_info& bi){

    assert(client_opt.is_bacon_allele || client_opt.is_bacon_snp);

    position_snp_call_bacon(pi,bi.bas);

    if(! bi.bas.is_valid) return;

    bi.unused_count = unused_count;

    const unsigned n_calls(pi.calls.size());
    for(unsigned i(0);i<n_calls;++i){
        assert(pi.calls[i].base_id!=BASE_ID::ANY);

        bi.base_count[pi.calls[i].base_id]++;
        bi.used_count++;
    }

    bi.is_max2_reported = (bi.bas.max2_score >= bacon_max2_score_print_thresh);

    // Calculate snp information even if only the bacon allele-caller
    // is requested to maintain compatibility with CASAVA's
    // "snp-calling by filtering the allele table" system. This way
    // anomaly filters will still be run on snp positions in the
    // allele table:
    //
    const unsigned ref_base_id(base_to_id(pi.ref_base));
    if(bi.bas.max_score >= client_opt.bacon_call_thresh) {
        const bool is_max_nonref(bi.bas.max_base_id != ref_base_id);
        const bool is_max2_sig(bi.bas.max2_score >= client_opt.bacon_second_call_thresh);
        bool is_max2_ratio_valid(false);
        if(is_max2_sig){
            is_max2_ratio_valid = ((bi.bas.max_score/bi.bas.max2_score) <= client_opt.bacon_het_snp_ratio_thresh);
        }
        if(is_max_nonref || (is_max2_sig && is_max2_ratio_valid)){
            bi.is_snp=true;
            bi.is_het=(is_max2_sig && is_max2_ratio_valid);
            bi.bst = bacon_snp_class(bi.bas.max_base_id,bi.bas.max2_base_id,ref_base_id,bi.is_het);
        }
    }
}



void
report_bacon_allele_call(const blt_streams& client_io,
                         const unsigned pos,
                         const bacon_info& bi){

    assert(client_io.bacon_allele_osptr());
    std::ostream& os(*client_io.bacon_allele_osptr());

    if(! bi.bas.is_valid) {
        static const char* null_allele_line = "\t0\t0\t0\t0\tX\t0\t0\t0.0\n";
        os << pos << null_allele_line;
        return;
    }


    os << pos << '\t';
    for(unsigned i(0);i<N_BASE;++i){ os << bi.base_count[i] << '\t'; }

    char max_base(id_to_base(bi.bas.max_base_id));
    if(bi.bas.max_score <= 0.) max_base = 'X';

    os << max_base;
    if(bi.is_max2_reported) os << id_to_base(bi.bas.max2_base_id);

    os << '\t'
       << bi.used_count+bi.unused_count << '\t'
       << bi.used_count << '\t';

    os.setf(std::ios::fixed,std::ios::floatfield);
    os << std::setprecision(2);
    os << bi.bas.max_score;
    if(bi.is_max2_reported) os << ':' << bi.bas.max2_score;
    os.unsetf(std::ios::fixed);

    os << '\n';
}



void
report_bacon_snp_call(const blt_streams& client_io,
                      const unsigned pos,
                      const char ref_base,
                      const bacon_info& bi){

    assert(client_io.bacon_snp_osptr() && bi.is_snp);

    std::ostream& os(*client_io.bacon_snp_osptr());

    os << pos << '\t';
    for(unsigned i(0);i<N_BASE;++i){ os << bi.base_count[i] << '\t'; }

    os << id_to_base(bi.bas.max_base_id);
    if(bi.is_het) os << id_to_base(bi.bas.max2_base_id);

    os << '\t'
       << bi.used_count+bi.unused_count << '\t'
       << bi.used_count << '\t';

    os.setf(std::ios::fixed,std::ios::floatfield);
    os << std::setprecision(2);
    os << bi.bas.max_score;
    if(bi.is_max2_reported) os << ':' << bi.bas.max2_score;
    os.unsetf(std::ios::fixed);

    os << '\t' << ref_base << '\t' << bacon_snp_label(bi.bst) << '\n';
}

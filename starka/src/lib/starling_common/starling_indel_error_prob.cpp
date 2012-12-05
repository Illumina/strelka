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


#include "blt_util/blt_exception.hh"
#include "blt_util/log.hh"
#include "starling_common/starling_indel_error_prob.hh"

#include <cmath>

#include <cstdlib>
#include <iostream>
#include <utility>



static
void
get_indel_error_prob_hpol_len(const unsigned hpol_len,
                              double& insert_error_prob,
                              double& delete_error_prob){

    // indel error model parameters for P(error) = Ax+Bx^C, where x=hpol_len
    // note that fit does not cover length 1 deletions,
    // for which the estimated value is instead provided directly
    //
    // \todo get these parameters out of the code!
    //
    static const double insert_A(5.03824e-7);
    static const double insert_B(3.30572e-10);
    static const double insert_C(6.99777);

    static const double delete_hpol1_err(3.00057e-6);
    static const double delete_A(1.09814e-5);
    static const double delete_B(5.19742e-10);
    static const double delete_C(6.99256);

    const double insert_g(insert_A*hpol_len+insert_B*std::pow(hpol_len,insert_C));
    insert_error_prob=(1.-std::exp(-insert_g));

    double delete_g(delete_hpol1_err);
    if(hpol_len>1){
        delete_g = delete_A*hpol_len+delete_B*std::pow(hpol_len,delete_C);
    }
    delete_error_prob=(1.-std::exp(-delete_g));
}



//
// "indel_error" is the probability that the read supporting the indel case is an error
// "ref_error" is the probability that the read supporting the ref case is an error
//
void
get_indel_error_prob(const starling_options& client_opt,
                     const starling_indel_report_info& iri,
                     double& indel_error_prob,
                     double& ref_error_prob){

    // cache results for any realistic homopolymer length:
    static const unsigned max_hpol_len(40);
    static bool is_init(false);

    // the pair is the spurious value for (insertion,deletion):
    //
    static std::pair<double,double> indel_error_prob_len[max_hpol_len];

    if(! is_init) {
        if(! client_opt.is_simple_indel_error) {
            double itmp(0);
            double dtmp(0);
            for(unsigned i(0);i<max_hpol_len;++i){
                get_indel_error_prob_hpol_len(i+1,itmp,dtmp);
                indel_error_prob_len[i] = std::make_pair(itmp,dtmp);
            }
        } else {
            const double ie(client_opt.simple_indel_error);
            for(unsigned i(0);i<max_hpol_len;++i){
                indel_error_prob_len[i] = std::make_pair(ie,ie);
            }
        }
        is_init=true;
    }

    const bool is_simple_indel(iri.it==INDEL::INSERT || iri.it==INDEL::DELETE);

    if(! is_simple_indel) {
        // breakpoints and swaps --
        // use zero repeat error for now.
        //
        // TODO - provide estimates for complex indels
        //
        indel_error_prob=std::max(indel_error_prob_len[0].first,indel_error_prob_len[0].second);
        ref_error_prob=indel_error_prob;
    } else {
        // treat everything besides simple homopolymer
        // contractions/expansions as homopolymer length 1:
        //
        if(iri.repeat_unit.size() == 1) {
            static const unsigned one(1);
            const unsigned ref_hpol_len = std::min(std::max(iri.ref_repeat_count,one),max_hpol_len);
            const unsigned indel_hpol_len = std::min(std::max(iri.indel_repeat_count,one),max_hpol_len);
            int indel_size(1);
            static const bool is_indel_size_dependent_error(false);
            if(is_indel_size_dependent_error) {
                indel_size=(std::abs(static_cast<long>(iri.ref_repeat_count)-
                                     static_cast<long>(iri.indel_repeat_count)));
            }

            if       (iri.it == INDEL::INSERT) {
                indel_error_prob=std::max(indel_error_prob_len[0].first,
                                          std::pow(indel_error_prob_len[ref_hpol_len-1].first,indel_size));
                ref_error_prob=std::max(indel_error_prob_len[0].second,
                                        std::pow(indel_error_prob_len[indel_hpol_len-1].second,indel_size));
            } else if(iri.it == INDEL::DELETE) {
                indel_error_prob=std::max(indel_error_prob_len[0].second,
                                          std::pow(indel_error_prob_len[ref_hpol_len-1].second,indel_size));
                ref_error_prob=std::max(indel_error_prob_len[0].first,
                                        std::pow(indel_error_prob_len[indel_hpol_len-1].first,indel_size));
            } else {
                log_os << "ERROR: Unknown indel type: " << iri.desc << "\n";
                throw blt_exception("Unknown indel type.");
            }
        } else {
            if       (iri.it == INDEL::INSERT) {
                indel_error_prob=indel_error_prob_len[0].first;
                ref_error_prob=indel_error_prob_len[0].second;
            } else if(iri.it == INDEL::DELETE) {
                indel_error_prob=indel_error_prob_len[0].second;
                ref_error_prob=indel_error_prob_len[0].first;
            } else {
                log_os << "ERROR: Unknown indel type: " << iri.desc << "\n";
                throw blt_exception("Unknown indel type.");
            }
        }
    }
}



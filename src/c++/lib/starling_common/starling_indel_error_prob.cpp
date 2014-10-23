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
                              double& delete_error_prob, const std::string& context)
{
    // Calculate p(error) of
    //    CASE: del
    //    FIT pars: [  1.49133831e-03   1.03348683e+01   1.13646811e+00   1.18488282e-05]
    //    Function prob(error)=0.00149133830825/ (1 + exp((10.3348683003-x)/1.13646810558))+1.18488281756e-05
    //    --------------------
    //    CASE: ins
    //    FIT pars: [  1.09573511e-03   9.82226042e+00   1.03579658e+00   8.31843836e-06]
    //    Function prob(error)=0.00109573511176/ (1 + exp((9.82226041538-x)/1.03579658224))+8.31843836296e-06
    //    --------------------

    // logistic model, fit for AT Reference context
    if (context=="AT")
    {
        static const double insert_A(1.49133831e-03);
        static const double insert_B(1.03348683e+01);
        static const double insert_C(1.13646811e+00);
        static const double insert_D(1.18488282e-05);

        static const double delete_A(1.09573511e-03);
        static const double delete_B(9.82226042e+00);
        static const double delete_C(1.03579658e+00);
        static const double delete_D(8.31843836e-06);

        const double insert_g(insert_A/ (1 + std::exp((insert_B-hpol_len)/insert_C))+insert_D);
        insert_error_prob=(1.-std::exp(-insert_g/hpol_len));

        const double delete_g(delete_A/ (1 + std::exp((delete_B-hpol_len)/delete_C))+delete_D);
        delete_error_prob=(1.-std::exp(-delete_g/hpol_len));
    }
    else
    {
        // new polynomial model, fit for CG reference context
        static const double insert_A(5.03824e-7);
        static const double insert_B(3.30572e-10);
        static const double insert_C(6.99777);

        static const double delete_hpol1_err(5.00057e-5);
        static const double delete_A(1.09814e-5);
        static const double delete_B(5.19742e-10);
        static const double delete_C(6.99256);

        const double insert_g(insert_A*hpol_len+insert_B*std::pow(hpol_len,insert_C));
        insert_error_prob=(1.-std::exp(-insert_g));

        double delete_g(delete_hpol1_err);
        if (hpol_len>1)
        {
            delete_g = delete_A*hpol_len+delete_B*std::pow(hpol_len,delete_C);
        }
        delete_error_prob=(1.-std::exp(-delete_g));
    }
}

// return the pre-calculated indel error rate for a given repeat-context and hpol length
static const unsigned max_hpol_len(40);
typedef std::pair<double,double> error_model[max_hpol_len];

struct PatternErrorModel
{
    PatternErrorModel()
    {
        const std::string AT_case = "AT";
        const std::string CG_case = "CG";

        double itmp(0);
        double dtmp(0);
        for (unsigned i(0); i<max_hpol_len; ++i)
        {
            get_indel_error_prob_hpol_len(i+1,itmp,dtmp,CG_case);
            indel_error_prob_len_CG[i] = std::make_pair(itmp,dtmp);
//            log_os << i << "_CG: " << itmp <<  " " << dtmp <<  "\n"; //print out test
            get_indel_error_prob_hpol_len(i+1,itmp,dtmp,AT_case);
            indel_error_prob_len_AT[i] = std::make_pair(itmp,dtmp);
//            log_os << i << "_AT: " << itmp <<  " " << dtmp <<  "\n"; //print out test
        }
    }

    const error_model&
    getModel(
        const std::string& overall_error_model,
        const std::string& pattern) const
    {
        // choose the error model based on
        if (overall_error_model=="old")
        {
            //        log_os << "Using indel error model: " << overall_error_model << "\n";
            return indel_error_prob_len_CG; // for now this is the old polynomial model
        }
        else if (overall_error_model=="stratified")
        {

            if ("G"==pattern || "C"==pattern)
            {
                return indel_error_prob_len_CG;
            }
            else
            {
                return indel_error_prob_len_AT;
            }
        }
        else
        {
            return indel_error_prob_len_AT;
        }
    }

    error_model indel_error_prob_len_AT;
    error_model indel_error_prob_len_CG;
};


static PatternErrorModel emodel;

#if 0
error_model&
get_pattern_error_model(
    const std::string& overall_error_model,
    std::string pattern="A",
    const int indel_length=1)
{

    // cache results for any realistic homopolymer length:
    // Treat everything above indel length 50 the same.
    // static const unsigned max_indel_len(50);
    static bool is_init(false);
    std::map<std::string,error_model> errorMap;
    // the pair is the spurious value for (insertion,deletion):
    // stratify by AT and CG case
    static error_model indel_error_prob_len_AT;
    static error_model indel_error_prob_len_CG;

    if (indel_length<0)
    {
        //TODO use indel length
    }

    // initialize error
    if (! is_init)
    {
        double itmp(0);
        double dtmp(0);
        for (unsigned i(0); i<max_hpol_len; ++i)
        {
            const std::string AT_case = "AT";
            const std::string CG_case = "CG";

            get_indel_error_prob_hpol_len(i+1,itmp,dtmp,CG_case);
            indel_error_prob_len_CG[i] = std::make_pair(itmp,dtmp);
//            log_os << i << "_CG: " << itmp <<  " " << dtmp <<  "\n"; //print out test
            get_indel_error_prob_hpol_len(i+1,itmp,dtmp,AT_case);
            indel_error_prob_len_AT[i] = std::make_pair(itmp,dtmp);
//            log_os << i << "_AT: " << itmp <<  " " << dtmp <<  "\n"; //print out test
        }
        is_init=true;
    }

    // choose the error model based on
    if (overall_error_model=="old")
    {
//        log_os << "Using indel error model: " << overall_error_model << "\n";
        return indel_error_prob_len_CG; // for now this is the old polynomial model
    }
    else if (overall_error_model=="stratified")
    {

        if ("G"==pattern or "C"==pattern)
        {
            return indel_error_prob_len_CG;
        }
        else
        {
            return indel_error_prob_len_AT;
        }
    }
    else
    {
        return indel_error_prob_len_AT;
    }
}
#endif


// "indel_error" is the probability that the read supporting the indel case is an error
// "ref_error" is the probability that the read supporting the ref case is an error
//
void
get_indel_error_prob(const starling_options& client_opt,
                     const starling_indel_report_info& iri,
                     double& indel_error_prob,
                     double& ref_error_prob)
{
    const bool is_simple_indel(iri.it==INDEL::INSERT || iri.it==INDEL::DELETE);

//    log_os << "Indel model " << client_opt.indel_error_model << "\n";

    const error_model& indel_error_prob_len(emodel.getModel(client_opt.indel_error_model,iri.repeat_unit));

    if (! is_simple_indel)
    {
        // breakpoints and swaps --
        // use zero repeat error for now.
        //
        // TODO - provide estimates for complex indels
        //
        indel_error_prob=std::max(indel_error_prob_len[0].first,indel_error_prob_len[0].second);
        ref_error_prob=indel_error_prob;
    }
    else
    {
        // treat everything besides simple homopolymer
        // contractions/expansions as homopolymer length 1:
        //
//        log_os << "Im doing in indels \n";
        if (iri.repeat_unit.size() == 1)
        {
            static const unsigned one(1);
            const unsigned ref_hpol_len = std::min(std::max(iri.ref_repeat_count,one),max_hpol_len);
            const unsigned indel_hpol_len = std::min(std::max(iri.indel_repeat_count,one),max_hpol_len);
            int indel_size(1);
            static const bool is_indel_size_dependent_error(false);
            if (is_indel_size_dependent_error)
            {
                indel_size=(std::abs(static_cast<long>(iri.ref_repeat_count)-
                                     static_cast<long>(iri.indel_repeat_count)));
            }

            if       (iri.it == INDEL::INSERT)
            {
                indel_error_prob=std::max(indel_error_prob_len[0].first,
                                          std::pow(indel_error_prob_len[ref_hpol_len-1].first,indel_size));
//            log_os << "error prob: " << indel_error_prob_len[ref_hpol_len-1].first << "\n";
                //reverse prob that true allele has been masked as reference by chance
                //may want to leave this term for now.
                ref_error_prob=client_opt.indel_ref_error_factor
                               * std::max(indel_error_prob_len[0].second,
                                        std::pow(indel_error_prob_len[indel_hpol_len-1].second,indel_size));
            }
            else if (iri.it == INDEL::DELETE)
            {
                indel_error_prob=std::max(indel_error_prob_len[0].second,
                                          std::pow(indel_error_prob_len[ref_hpol_len-1].second,indel_size));
                ref_error_prob=client_opt.indel_ref_error_factor
                               * std::max(indel_error_prob_len[0].first,
                                        std::pow(indel_error_prob_len[indel_hpol_len-1].first,indel_size));
            }
            else
            {
                log_os << "ERROR: Unknown indel type: " << iri.desc << "\n";
                throw blt_exception("Unknown indel type.");
            }
        }
        else
        {
            if (iri.it == INDEL::INSERT)
            {
                indel_error_prob=indel_error_prob_len[0].first;
                ref_error_prob=indel_error_prob_len[0].second;
            }
            else if (iri.it == INDEL::DELETE)
            {
                indel_error_prob=indel_error_prob_len[0].second;
                ref_error_prob=indel_error_prob_len[0].first;
            }
            else
            {
                log_os << "ERROR: Unknown indel type: " << iri.desc << "\n";
                throw blt_exception("Unknown indel type.");
            }
        }
    }
}



/*
 * cmodel.cpp
 *
 *  Created on: Dec 1, 2013
 *      Author: Morten Kallberg
 */

#include "cmodel.hh"
#include <stdlib.h>     /* atof */
#include <iostream>
#include <sstream>
#include <cassert>
#include <string>
#include <fstream>
#include <iterator>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include "blt_util/qscore.hh"

//#define DEBUG_MODEL

#ifdef DEBUG_MODEL
#include "blt_util/log.hh"
#endif


// add model paramaters
void c_model::add_parameters(const parmap& myPars) {
    this->pars = myPars;
}

void c_model::do_rule_model(featuremap& cutoffs, site_info& si) {
    if (si.smod.gqx<cutoffs["GQX"]) si.smod.set_filter(VCF_FILTERS::LowGQX);
    if (cutoffs["DP"]>0) {
        if ((si.n_used_calls+si.n_unused_calls) > cutoffs["DP"]) si.smod.set_filter(VCF_FILTERS::HighDepth);
    }

    // high DPFratio filter
    const unsigned total_calls(si.n_used_calls+si.n_unused_calls);
    if (total_calls>0) {
        const double filt(static_cast<double>(si.n_unused_calls)/static_cast<double>(total_calls));
        if (filt>cutoffs["DPFratio"]) si.smod.set_filter(VCF_FILTERS::HighBaseFilt);
    }

    if (si.dgt.is_snp) {
        if (si.dgt.sb>cutoffs["HighSNVSB"]) si.smod.set_filter(VCF_FILTERS::HighSNVSB);
    }
}

//Transform the features with the specified scaling parameters that were used to standardize
//the dataset to zero mean and unit variance: newVal = (oldVal-centerVal)/scaleVal.
featuremap c_model::normalize(featuremap features, featuremap& adjust_factor, featuremap& norm_factor) {
    for (featuremap::const_iterator it = norm_factor.begin(); it != norm_factor.end(); ++it) { // only normalize the features that are needed
        features[it->first] = (features[it->first]-adjust_factor[it->first])/norm_factor[it->first];
    }
    return features;
}

//Model: ln(p(TP|x)/p(FP|x))=w_1*x_1 ... w_n*x_n + w_1*w_2*x_1 ... w_{n-1}*w_{n}*x_1
// calculate sum from feature map
double c_model::log_odds(featuremap features, featuremap& coeffs) {
    using namespace boost::algorithm;
    std::vector<std::string> tokens;
    double sum = coeffs["Intercept"];
//    log_os << "sum" << "=" << sum << "\n";
    for (featuremap::const_iterator it = coeffs.begin(); it != coeffs.end(); ++it) {
        if (it->first !="Intercept" && it->second !=0) { // check that our coefficient is different from 0
            split(tokens, it->first, is_any_of(":"));
//            log_os << it->first << "=" << it->second << "\n";
            double term = it->second;
            for (unsigned int i=0; i < tokens.size(); i++) {
                term = term*features[tokens[i]];
//                log_os << tokens[i] << "=" << features[tokens[i]] << "\n";
            }
//            log_os << "term" << "=" << term << "\n";
//            log_os << "\n";
            // use term to determine the most predictive parameter
            sum += term;
//            log_os << "sum " << "=" << sum << "\n";
//            log_os << tokens.size() << "\n";
        }
    }
    //TO-DO sort map to find most predictive feature for setting filter for
    return sum;
}

//- From the identities logOddsRatio = ln(p(TP|x)/p(FP|x)) and p(TP|x) + p(FP|x) = 1,
//  solve for p(TP): p(TP) = 1/(1+exp(-logOddsRatio))
//- Rescale the probabilities p(TP), p(FP) with the specified class priors,
//  as they were counted on the training set before balancing to 50/50.
//- Rescale with the class priors; hard-coded 0.5, because we have 1:1 ratio in the balanced data
//  p(FP).pos.res = p(FP).pos*(minorityPrior/0.5)
//  p(TP).pos.res = p(TP).pos*(minorityPrior/0.5)
//- Renormalize to sum to one
//  p(FP).pos.ret = p(FP).pos.res/(p(TP).pos.res + p(FP).pos.res)
//  p(TP).pos.ret = p(TP).pos.res/(p(TP).pos.res + p(FP).pos.res)
//- Above simplifies to:
//  p(FP).pos.ret = p(FP).pos*minorityPrior/(1+2*minorityPrior*p(FP).pos-minorityPrior-p(FP).pos)
//- Convert the rescaled probability p(FP) into a
//  Q-score: Q-score = round( -10*log10(p(FP).pos.ret) )
static
int prior_adjustment(
    const double raw_score,
    const double minorityPrior) {

    double pFP          = 1.0/(1+std::exp(raw_score)); // this calculation can likely be simplified
    double pFPrescale   = pFP*minorityPrior/(1+2*minorityPrior*pFP-minorityPrior-pFP);
    int qscore          = error_prob_to_qphred(pFPrescale);
    #ifdef DEBUG_MODEL
        log_os << "minorityPrior " << minorityPrior << "\n";
        log_os << "pFP=" << pFP << "\n";
        log_os << "rescale=" << pFPrescale << "\n";
    //        log_os << "experimental=" << qscore_test << "\n";
    #endif

    // cap the score at 40
    if (qscore>40)
        qscore = 40;
    // TODO check for inf and NaN artifacts

    return qscore;
}
void c_model::apply_qscore_filters(site_info& si, featuremap& qscore_cuts){//, featuremap& most_predictive) {
//    most_predictive.size();
    if (si.Qscore < qscore_cuts["Q"]) {
        si.smod.set_filter(VCF_FILTERS::LowGQX); // more sophisticated filter setting here
    }
}

void c_model::score_instance(featuremap features, site_info& si) {
    if (this->model_type=="LOGISTIC") { //case we are using a logistic regression mode
        std::string snpCase = "homsnp";
        if (si.is_het())
            snpCase = "hetsnp";
//        this->sanity_check();
        // normalize
        featuremap norm_features = this->normalize(features,this->pars[snpCase]["scalecenter"],this->pars[snpCase]["scaleshift"]);

        //calculates log-odds ratio
        double raw_score = this->log_odds(norm_features,this->pars[snpCase]["coefs"]);

        // adjust by prior and calculate q-score
        si.Qscore = prior_adjustment(raw_score,this->pars[snpCase]["priors"]["minorityPrior"]);

        // set filters according to q-scores
//        featuremap most_pred; //place-holder
        this->apply_qscore_filters(si,this->pars[snpCase]["qcutoff"]);

#ifdef DEBUG_MODEL
//        log_os << "Im doing a logistic model" << "\n";
#endif

    }
    else if (this->model_type=="RULE") { //case we are using a rule based model
        this->do_rule_model(this->pars["snp"]["cutoff"],si);
    }
}

// TODO decompose to unit-test
//void c_model::sanity_check(){
//    featuremap fm;
//    featuremap fm2;
//    fm["GQX"]   = 42;
//    fm["DP"]    = 28;
//    fm["AD2"]   = 6;
//    fm["SNVSB"] = -5.9;
//    fm["SNVHPOL"] = 3;
//    fm["VFStar"] = 0.214286;
//    fm["DPF"] = 0;
//    fm["MQ"] = 60;
//
//    fm2["GQX"]   = 128;
//    fm2["DP"]    = 34;
//    fm2["AD2"]   = 12;
//    fm2["SNVSB"] = -18;
//    fm2["SNVHPOL"] = 3;
//    fm2["VFStar"] = 0.333333;
//    fm2["DPF"] = 2;
//    fm2["MQ"] = 60;
//
//    std::string snpCase = "hetsnp";
//    featuremap norm_features = this->normalize(fm2,this->pars[snpCase]["scalecenter"],this->pars[snpCase]["scaleshift"]);
////    for (featuremap::const_iterator it = norm_features.begin(); it != norm_features.end(); ++it) {
////        log_os << it->first << "=" << it->second << "\n";
////    }
//    const double raw_score = this->log_odds(norm_features,this->pars[snpCase]["coefs"]);
//    log_os << "score " << raw_score << std::endl;
//    int q = prior_adjustment(raw_score,this->pars[snpCase]["priors"]["minorityPrior"]);
//    log_os << "q " << q << std::endl;
//
//}


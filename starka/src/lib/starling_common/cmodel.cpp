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
#include <exception>
#include <string>
#include <fstream>
#include <iterator>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <stdlib.h>     /* atof */

#define DEBUG_MODEL

#ifdef DEBUG_MODEL
#include "blt_util/log.hh"
#endif


c_model::c_model(std::string name, std::string type){
    this->model_name = name;
    this->model_type = type;
}
c_model::~c_model() {
    // TODO Auto-generated destructor stub
}

// add model paramaters
void c_model::add_parameters(parmap myPars){
    this->pars = myPars;
}


void c_model::do_rule_model(featuremap& cutoffs, site_info& si){
      if (si.smod.gqx<cutoffs["GQX"]) si.smod.set_filter(VCF_FILTERS::LowGQX);
//      log_os << "GQX set to " << cutoffs["DPFratio"] << "\n";
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
featuremap c_model::normalize(featuremap features, featuremap& adjust_factor, featuremap& norm_factor){
    for(featuremap::const_iterator it = norm_factor.begin(); it != norm_factor.end(); ++it){ // only normalize the features that are needed
//        log_os << it->first << "=" << it->second << " ";
//        log_os << "scale " << "=" << adjust_factor[it->first] << " ";
        features[it->first] = (features[it->first]-adjust_factor[it->first])/norm_factor[it->first];
//        log_os << it->first << "=" << it->second << "\n";
    }
    return features;
}

//Model: ln(p(TP|x)/p(FP|x))=w_1*x_1 ... w_n*x_n + w_1*w_2*x_1 ... w_{n-1}*w_{n}*x_1
// calculate sum from feature map
double c_model::log_odds(featuremap features, featuremap& coeffs){
    using namespace boost::algorithm;
    std::vector<std::string> tokens;
    std::map<std::string,double> predictive;
    double sum = coeffs["Intercept"];
//    log_os << "sum" << "=" << sum << "\n";
    for(featuremap::const_iterator it = coeffs.begin(); it != coeffs.end(); ++it){
        if (it->first !="Intercept" && it->second !=0){ // check that out coefficient is greater than 0
            split(tokens, it->first, is_any_of(":"));
//            log_os << it->first << "=" << it->second << "\n";
            double term = it->second;
//            log_os << "term" << "=" << term << "\n";
            for (unsigned int i=0; i < tokens.size(); i++) {
                term = term*features[tokens[i]];
//                log_os << "term" << "=" << term << "\n";
            }
            // use term to determine the most predictive parameter
            predictive[it->first] = term;
            sum += term;
//            log_os << "sum " << "=" << sum << "\n";
//            log_os << tokens.size() << "\n";
        }
    }
    //TO-DO sort map to find most predictive feature for setting filter for
    return sum;
}

//- From the identities logOddsRatio = ln(p(TP|x)/p(FP|x)) and p(TP|x) + p(FP|x) = 1,
//solve for p(TP): p(TP) = 1/(1+exp(-logOddsRatio))
//- Rescale the probabilities p(TP), p(FP) with the specified class priors,
//as they were counted on the training set before balancing to 50/50.
//- Convert the rescaled probability p(FP) into a q score
//Q-score: Q-score = round( -10*log10(p(FP)) )
// simplification possible  qscore(raw) = round(10log((1+e^raw)/prior))
double c_model::prior_adjustment(const double raw_score, featuremap& priors){
    double pFP = 1.0 - 1.0/(1+exp(-raw_score)); // this calculation can likely be simplified
    double pFPrescale   = pFP*priors["minorityPrior"];
    double qscore       = round(-10.0 * log10(pFPrescale));
//    double qscore_test  = round(10*log10((1+exp(raw_score))/priors["minorityPrior"]));
    #ifdef DEBUG_MODEL
        log_os << "minorityPrior " << priors["minorityPrior"] << "\n";
        log_os << "pFP=" << pFP << "\n";
        log_os << "rescale=" << pFPrescale << "\n";
//        log_os << "experimental=" << qscore_test << "\n";
    #endif

    return qscore;
}



void c_model::score_instance(featuremap features, site_info& si){
    if (this->model_type=="LOGISTIC"){ //case we are using a logistic regression mode
        std::string snpCase = "homsnp";
        if (si.is_het())
            snpCase = "hetsnp";

        // normalize
        featuremap norm_features = this->normalize(features,this->pars[snpCase]["scalecenter"],this->pars[snpCase]["scaleshift"]);

        //calculates log-odds ratio
        double raw_score = this->log_odds(norm_features,this->pars[snpCase]["coefs"]);

        // adjust by prior and calculate q-score
        si.Qscore = this->prior_adjustment(raw_score,this->pars[snpCase]["priors"]);
    #ifdef DEBUG_MODEL
            log_os << "Im doing a logistic model" << "\n";
    #endif

    }
    else if (this->model_type=="RULE"){//case we are using a rule based model
        featuremap myCutoffs;
        this->do_rule_model(this->pars["snp"]["cutoff"],si);
    }
    else{                             //should never end up here
        si.Qscore = 10.1;
    }
}


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

//#define DEBUG_MODEL

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
void c_model::add_parameter(std::vector<std::string> tokens,std::string submodel, std::string context){
    if (tokens.size()>1){
//        log_os << "Loaded " << submodel << " " << context << " " << tokens.at(0) << " " << tokens.at(1) << "\n";
        this->pars[submodel][context].insert(std::make_pair(tokens.at(0),atof(tokens.at(1).c_str())));
//        log_os << this->pars[submodel][context][tokens.at(0)] << "\n";
//
//        log_os << "map size " << this->pars[submodel][context].size() << "\n";
    }
}

void c_model::do_rule_model(featuremap cutoffs, site_info& si){
      if (si.smod.gqx<cutoffs["GQX"]) si.smod.set_filter(VCF_FILTERS::LowGQX);
//      log_os << "GQX set to " << cutoffs["GQX"] << "\n";
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

          //kick out hpol filter completely
//          if (opt.is_max_snv_hpol) {
//              if (static_cast<int>(si.hpol)>opt.max_snv_hpol) si.smod.set_filter(VCF_FILTERS::HighSNVHPOL);
//          }
      }
}


featuremap c_model::normalize(featuremap features, featuremap adjust_factor, featuremap norm_factor){
    return features;
}

double c_model::log_odds(featuremap features, featuremap scaling_factor){
    return 0.08934589;
}

// adjust for prior and calculate q-score
double c_model::prior_adjustment(double raw_score){
    return 0.18934589;
}



void c_model::score_instance(featuremap features, site_info& si){
    if (this->model_type=="LOGISTIC"){ //case we are using a logistic regression mode
        // pull out the correct model parameters for the case we are dealing with

        // normalize
        featuremap norm_features = this->normalize(features,features,features);

        //calculates log-odds ratio
        double raw_score = this->log_odds(norm_features,features);

        // adjust by prior and calculate q-score
        si.Qscore = this->prior_adjustment(raw_score);
    #ifdef DEBUG_MODEL
            log_os << "Im doing a logistic model" << "\n";
    #endif

    }
    else if (this->model_type=="RULE"){//case we are using a rule based model
        featuremap myCutoffs;
        myCutoffs["GQX"] = 30;
        myCutoffs["HighSNVSB"] = 10;
        myCutoffs["DPFratio"]  = 0.6;
        myCutoffs["DP"]  = 0;
        this->do_rule_model(myCutoffs,si);
    }
    else{                             //should never end up here
        si.Qscore = 10.1;
    }
}


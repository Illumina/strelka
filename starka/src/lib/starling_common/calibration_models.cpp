/*
 * calbrationmodels.cpp
 *
*  Created on: Oct 10, 2013
 * Author: Morten Kallberg
 */

#include "calibration_models.hh"
//#include <boost/property_tree/ptree.hpp>
//#include <boost/property_tree/json_parser.hpp>
//#include <boost/foreach.hpp>
#include <cassert>
#include <exception>
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <iterator>
#include <map>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>


//#define DEBUG_CAL

#ifdef DEBUG_CAL
#include "blt_util/log.hh"
#endif

calibration_models::calibration_models() {
    this->set_model("default");
}
calibration_models::~calibration_models(){};

void calibration_models::clasify_site(const gvcf_options& opt, const gvcf_deriv_options& dopt, site_info& si){
    // create site value feature dict
    featuremap features = si.get_qscore_features();

    if (si.dgt.is_snp && this->model_name!="default") {
//        for(featuremap::const_iterator it = features.begin(); it != features.end(); ++it)
//            log_os << it->first << "=" << it->second << " ";
//        log_os << "\n";

        c_model myModel = this->get_model(this->model_name);
        myModel.score_instance(features,si);
    }
    else{
        // don't know what to do with this site, throw it to the old default filters
        this->default_clasify_site(opt,dopt,si);
    }
}

void calibration_models::default_clasify_site(const gvcf_options& opt, const gvcf_deriv_options& dopt, site_info& si){
    if (opt.is_min_gqx) {
          if (si.smod.gqx<opt.min_gqx) si.smod.set_filter(VCF_FILTERS::LowGQX);
      }

      if (dopt.is_max_depth) {
          if ((si.n_used_calls+si.n_unused_calls) > dopt.max_depth) si.smod.set_filter(VCF_FILTERS::HighDepth);
      }

      // high DPFratio filter
      if (opt.is_max_base_filt) {
          const unsigned total_calls(si.n_used_calls+si.n_unused_calls);
          if (total_calls>0) {
              const double filt(static_cast<double>(si.n_unused_calls)/static_cast<double>(total_calls));
              if (filt>opt.max_base_filt) si.smod.set_filter(VCF_FILTERS::HighBaseFilt);
          }
      }

      if (si.dgt.is_snp) {
          if (opt.is_max_snv_sb) {
              if (si.dgt.sb>opt.max_snv_sb) si.smod.set_filter(VCF_FILTERS::HighSNVSB);
          }

          if (opt.is_max_snv_hpol) {
              if (static_cast<int>(si.hpol)>opt.max_snv_hpol) si.smod.set_filter(VCF_FILTERS::HighSNVHPOL);
          }
      }
}

void calibration_models::set_model(const std::string& name){
    modelmap::iterator it = this->models.find(name);
    if(it != this->models.end())
    {
        this->model_name = name;
    }
    else{
        this->model_name = "default";
    }
    #ifdef DEBUG_CAL
        log_os << "Calibration model set to '" << this->model_name << "'\n";
    #endif
}

c_model calibration_models::get_model(std::string name){
    modelmap::iterator it;
    it = this->models.find(name);
    return it->second;
}


void calibration_models::load_models(std::string model_file){
    using namespace boost::algorithm;

    #ifdef DEBUG_CAL
        log_os << "Loading models from file: " << model_file << "\n";
    #endif
    std::ifstream myReadFile;
    myReadFile.open(model_file.c_str());
    std::string output;
    std::string parspace;
    std::string submodel;
    std::string current_name;
    if (myReadFile.is_open()) {
        while (!myReadFile.eof()) {
           std::getline (myReadFile,output);
           std::vector<std::string> tokens;
           split(tokens, output, is_any_of(" ")); // tokenize string
           //case new model
           if (tokens.at(0).substr(0,3)=="###"){
               current_name = tokens.at(1);
               c_model current_model(tokens.at(1),tokens.at(2));
               this->models.insert(modelmap::value_type(current_name, current_model));
               #ifdef DEBUG_CAL
                   log_os << "Loading model: " << tokens.at(1) << " Type: " << tokens.at(2) << "\n";
               #endif
           }
           //load submodel
           else if(tokens.at(0)=="#"){
               #ifdef DEBUG_CAL
                   log_os << "submodel: " << tokens.at(1) << " parspace: " << tokens.at(2) << "\n";
               #endif
               submodel = tokens.at(1);
               parspace = tokens.at(2);
           }
           //case load parameters
           else{
               this->get_model(current_name).add_parameter(tokens,submodel,parspace);
           }
        }
    }
    #ifdef DEBUG_CAL
        log_os << "Done loading models" << "\n";
    #endif
    myReadFile.close();
}





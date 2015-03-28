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
/*
 * scoringmodels.cpp
 *
 *  Created on: Aug 20, 2014
 *      Author: mkallberg
 */
#include "scoringmodels.hh"

#include "blt_util/log.hh"
#include "blt_util/qscore.hh"
#include "blt_util/parse_util.hh"

#include "boost/property_tree/json_parser.hpp"

#include <cassert>

#include <iostream>
#include <sstream>

using boost::property_tree::ptree;
using namespace illumina::blt_util;

//#define DEBUG_SCORINGMODELS

#ifdef DEBUG_SCORINGMODELS
#include "blt_util/log.hh"
#endif



void indel_model::add_prop(const unsigned hpol_case, const double prop_ins,const double prop_del)
{
    if (hpol_case>0 && hpol_case<max_hpol_len)
    {
        this->model[hpol_case-1] = std::make_pair(prop_ins,prop_del);

    }
}

double indel_model::get_prop(const unsigned hpol_case) const
{
    return model[std::min(hpol_case,max_hpol_len)-1].first;
}



// Global static pointer used to ensure a single instance of the class.
scoring_models* scoring_models::m_pInstance = nullptr;

/** This function is called to create an instance of the class.
    Calling the constructor publicly is not allowed. The constructor
    is private and is only called by this Instance function.
*/
scoring_models& scoring_models::Instance()
{
    if (!m_pInstance)   // Only allow one instance of class to be generated.
        m_pInstance = new scoring_models;
    return *m_pInstance;
}


double scoring_models::score_instance(const feature_type& features) const
{
    const double score = this->randomforest_model.getProb(features);
    return error_prob_to_phred(score);
}

const error_model& scoring_models::get_indel_model(const std::string& /*pattern*/) const
{
    return this->indel_models.at(this->current_indel_model).model;
}


void scoring_models::load_indel_model(const ptree& pt,const std::string& model_name)
{

    std::string s = imodels + "." + model_name;
//    log_os << s << std::endl;
    indel_model temp_model;
    unsigned i=0;
    for (const ptree::value_type& v : pt.get_child(s))
    {
        temp_model.add_prop(i,atof(v.second.data().c_str()),atof(v.second.data().c_str()));
        i++;
    }
    this->indel_models[model_name] = temp_model;
    this->indel_init = true;
}



void scoring_models::load_calibration_model(const ptree& pt,const std::string& model_name,const std::string& model_type)
{
    if (model_name!="" && model_type!="")
    {
//        log_os << "Loading cali model: " << model_name << std::endl;
    }

    //TODO add case here for indels, currently only loading snps. Need to add another nesting level here
    this->randomforest_model.load(pt.get_child(model_name));
    this->calibration_init = true;
}

void scoring_models::load_models(const std::string& model_file)
{
    // assume file exists has been checked

    std::stringstream ss;
    std::ifstream file( model_file );
    ss << file.rdbuf();
    file.close();

    ptree pt;
    boost::property_tree::read_json(ss, pt);

    //load indel models
    for (const ptree::value_type& v : pt.get_child(imodels))
    {
//         log_os << "Reading indel model " << v.first <<  std::endl;
        this->load_indel_model(pt,v.first);
    }

    //load calibration models
    for (const ptree::value_type& v : pt.get_child(cmodels))
    {
        this->load_calibration_model(pt.get_child(cmodels),v.first);
    }
}

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
 * scoringmodels.hh
 *
 *  Created on: Aug 20, 2014
 *      Author: mkallberg
 */

#pragma once

#include "calibration/RandomForestModel.hh"

#include "boost/property_tree/ptree.hpp"

#include <map>
#include <string>
#include <vector>


static const std::string imodels="IndelModels";
static const std::string cmodels="CalibrationModels";
static const unsigned max_hpol_len(40);
static const unsigned max_indel_len(15);
typedef std::pair<double,double> error_model[max_hpol_len];

struct indel_model
{
    void add_prop(const unsigned hpol_case, const double prop_ins,const double prop_del);
    double get_prop(const unsigned hpol_case) const;
    error_model model;
private:
    std::string name;
};



struct scoring_models
{
    static scoring_models& Instance();
    void load_models(const std::string& model_file);
    void load_indel_model(const boost::property_tree::ptree& pt,const std::string& model_name);

    void load_calibration_model(const boost::property_tree::ptree& pt,const std::string& model_name,const std::string& model_type="RF");
    double score_instance(const feature_type& features) const;

    const error_model& get_indel_model(const std::string& pattern) const;
    bool indel_init=false;

    bool calibration_init=false;

private:
    scoring_models() {} // Private so that it can  not be called
    scoring_models(scoring_models const&) {}            // copy constructor is private
    scoring_models& operator=(scoring_models const&);  // assignment operator is private
    static scoring_models* m_pInstance;
    typedef std::map<std::string,indel_model> indel_modelmap;
    indel_modelmap indel_models;
    std::string current_indel_model;

    //typedef std::map<int,calibration_model> calibration_modelmap;
    //calibration_modelmap calibration_models;
    RandomForestModel randomforest_model;
};

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

#include "boost/property_tree/ptree.hpp"

#include <map>
#include <string>
#include <vector>



namespace STRELKA_VQSR_FEATURES
{

enum index_t
{
    QSS_NT,
    N_FDP_RATE,
    T_FDP_RATE,
    N_SDP_RATE,
    T_SDP_RATE,
    N_DP_RATE,
    TIER1_ALLELE_RATE,
    MQ,
    n_mapq0,
    strandBias,
    ReadPosRankSum,
    altmap,
    altpos,
    pnoise,
    pnoise2,
    SIZE
};

inline
const char*
get_feature_label(const unsigned idx)
{
    switch (idx)
    {
    case QSS_NT:
        return "QSS_NT";
    case N_FDP_RATE:
        return "N_FDP_RATE";
    case T_FDP_RATE:
        return "T_FDP_RATE";
    case N_SDP_RATE:
        return "N_SDP_RATE";
    case T_SDP_RATE:
        return "T_SDP_RATE";
    case N_DP_RATE:
        return "N_DP_RATE";
    case TIER1_ALLELE_RATE:
        return "TIER1_ALLELE_RATE";
    case MQ:
        return "MQ";
    case n_mapq0:
        return "n_mapq0";
    case strandBias:
        return "strandBias";
    case ReadPosRankSum:
        return "ReadPosRankSum";
    case altmap:
        return "altmap";
    case altpos:
        return "altpos";
    case pnoise:
        return "pnoise";
    case pnoise2:
        return "pnoise2";
    default:
        assert(false && "Unknown feature");
        return nullptr;
    }
}
};

typedef std::map<int, double> feature_type;


//namespace CALIBRATION_MODEL
//{
//
//    enum model_case
//    {
//        HetSNP,
//    };
//
//    inline
//    const char*
//    get_label(const unsigned idx)
//    {
//        switch (idx)
//        {
//        case HetSNP:
//            return "snphet";
//        default:
//            assert(0);
//            return NULL;
//        }
//    }
//}

static const std::string imodels="IndelModels";
static const std::string cmodels="CalibrationModels";
static const unsigned max_hpol_len(40);
static const unsigned max_indel_len(15);
typedef std::pair<double,double> error_model[max_hpol_len];

class indel_model
{
public:
    indel_model() {}
    void add_prop(const unsigned hpol_case, const double prop_ins,const double prop_del);
    double get_prop(const unsigned hpol_case) const;
    error_model model;
private:
    std::string name;
};



class calibration_model
{
public:
    calibration_model() {}

    typedef std::map<int, std::vector<double> > calibration_type;
    typedef std::vector< calibration_type > set_of_calibrations_type;

    void populate_storage_metadata();
    void load(const boost::property_tree::ptree& pt);

    double get_randomforest_proba(const feature_type& features) const;
private:
    double get_single_dectree_proba(const feature_type& features, int tree_index) const;

    int n_trees = 10;
    std::vector<std::string> calibration_data_names;
    set_of_calibrations_type all_trees;
    set_of_calibrations_type all_node_votes;
    set_of_calibrations_type all_decisions;
    std::vector< set_of_calibrations_type > all_rf_json_data;

    std::string name;
};

class scoring_models
{
public:
    static scoring_models* Instance();
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
    calibration_model randomforest_model;


    //
    //std::map< std::map<int, std::vector<double> >, double> proba_for_features;


};


/*
 * scoringmodels.hh
 *
 *  Created on: Aug 20, 2014
 *      Author: mkallberg
 */


#include <stdlib.h>     /* atof */
#include <iostream>
#include <sstream>
#include <cassert>
#include <string>
//#include <fstream>
//#include <iterator>
#include <map>
#include <vector>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>
#ifndef SCORINGMODELS_HH_
#define SCORINGMODELS_HH_

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
typedef std::map<int, double> feature_type;

class indel_model{
public:
    indel_model(){};
    void add_prop(const unsigned hpol_case, const double prop_ins,const double prop_del);
    double get_prop(const unsigned hpol_case);
    error_model model;
private:
    std::string name;
};

class calibration_model{
public:
    calibration_model(){};

    typedef std::map<int, std::vector<double> > calibration_type;
    typedef std::vector< calibration_type > set_of_calibrations_type;

    void populate_storage_metadata();
    void load(boost::property_tree::ptree pt);

    double get_randomforest_proba(const feature_type& features);
private:
    double get_single_dectree_proba(const feature_type& features, int tree_index);

    int n_trees = 10;
    std::vector<std::string> calibration_data_names;
    set_of_calibrations_type all_trees;
    set_of_calibrations_type all_node_votes;
    set_of_calibrations_type all_decisions;
    std::vector< set_of_calibrations_type > all_rf_json_data;



private:
    std::string name;
};

class scoring_models{
public:
   static scoring_models* Instance();
   void load_models(const std::string& model_file);
   void load_indel_model(boost::property_tree::ptree pt,const std::string& model_name);

   void load_calibration_model(boost::property_tree::ptree pt,const std::string& model_name,const std::string& model_type="RF");
   double score_instance(const feature_type& features);

   error_model& get_indel_model(const std::string& pattern);
   bool indel_init=false;

   bool calibration_init=false;

private:
   scoring_models(){};  // Private so that it can  not be called
   scoring_models(scoring_models const&){};             // copy constructor is private
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


#endif /* SCORINGMODELS_HH_ */

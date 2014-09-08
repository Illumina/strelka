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

static const unsigned max_hpol_len(40);
static const unsigned max_indel_len(15);

class indel_model{
public:
    indel_model(){};
    void add_prop(const unsigned hpol_case, const double prop);
    double get_prop(const unsigned hpol_case);
private:
    double emodel[max_hpol_len];
    std::string name;
};

class scoring_models{
public:
   static scoring_models* Instance();
   void load_models(const std::string& model_file);
   void load_indel_models(boost::property_tree::ptree pt,const std::string model_name);

private:
   scoring_models(){};  // Private so that it can  not be called
   scoring_models(scoring_models const&){};             // copy constructor is private
   scoring_models& operator=(scoring_models const&);  // assignment operator is private
   static scoring_models* m_pInstance;
   typedef std::map<std::string,indel_model> modelmap;
   modelmap models;
};


#endif /* SCORINGMODELS_HH_ */

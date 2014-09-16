/*
 * scoringmodels.cpp
 *
 *  Created on: Aug 20, 2014
 *      Author: mkallberg
 */
#include "blt_util/scoringmodels.hh"

#define DEBUG_SCORINGMODELS

#ifdef DEBUG_SCORINGMODELS
    #include "blt_util/log.hh"
#endif



// Global static pointer used to ensure a single instance of the class.
scoring_models* scoring_models::m_pInstance = nullptr;

/** This function is called to create an instance of the class.
    Calling the constructor publicly is not allowed. The constructor
    is private and is only called by this Instance function.
*/
scoring_models* scoring_models::Instance()
{
   if (!m_pInstance)   // Only allow one instance of class to be generated.
      m_pInstance = new scoring_models;
   return m_pInstance;
}

void indel_model::add_prop(const unsigned hpol_case, const double prop_ins,const double prop_del){
    if (hpol_case>0 && hpol_case<max_hpol_len){
        this->model[hpol_case-1] = std::make_pair(prop_ins,prop_del);
        #ifdef DEBUG_SCORINGMODELS
//            log_os << "Adding case " << hpol_case << " with prop " << prop << std::endl;
        #endif
    }
 }

double indel_model::get_prop(const unsigned hpol_case){
    if (hpol_case>40)
        return 1.0;
    return 1.0;
//        return this->model[max_hpol_len-1][0];
//    return this->model[hpol_case-1][0];

}

error_model& scoring_models::get_indel_model(const std::string& pattern){
    if (pattern=="f"){

    }
    return this->indel_models[this->current_indel_model].model;
}

double scoring_models::score_instance(const std::map<std::string,double> features){

    double score = 0.5;
    if(features.empty()){

    }

    return score; // returns calibration score
}


void scoring_models::load_indel_models(boost::property_tree::ptree pt,const std::string model_name){
    std::string s = imodels + "." + model_name;
//    log_os << s << std::endl;
    indel_model temp_model;
    unsigned i=0;
    BOOST_FOREACH(boost::property_tree::ptree::value_type &v, pt.get_child(s))
    {
        log_os << "Adding case " << i << " with prop " << v.second.data().c_str() << std::endl;
        temp_model.add_prop(i,atof(v.second.data().c_str()),atof(v.second.data().c_str()));
        i++;
    }
    this->indel_models[model_name] = temp_model;
    this->indel_init = true;
}

//void scoring_models::load_calibration_models(boost::property_tree::ptree pt,const std::string model_name){
//    if(model_name==""){}
//
//    this->calibration_init = true;
//}

void scoring_models::load_models(const std::string& model_file){
    // assume file exists has been checked
    std::stringstream ss;
    std::ifstream file( model_file );
    ss << file.rdbuf();
    file.close();

     boost::property_tree::ptree pt;
     boost::property_tree::read_json(ss, pt);

     //load indel models
     BOOST_FOREACH(boost::property_tree::ptree::value_type &v, pt.get_child(imodels))
         this->load_indel_models(pt,v.first);

     //load calibration models
//     BOOST_FOREACH(boost::property_tree::ptree::value_type &v, pt.get_child(cmodels))
//         this->load_calibration_models(pt,v.first);

}

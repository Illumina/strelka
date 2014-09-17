/*
 * scoringmodels.cpp
 *
 *  Created on: Aug 20, 2014
 *      Author: mkallberg
 */
#include "scoringmodels.hh"

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

    }
 }

double indel_model::get_prop(const unsigned hpol_case){
    if (hpol_case>40)
        return this->model[max_hpol_len-1].first;
    return this->model[hpol_case-1].first;
}

void calibration_model::populate_storage_metadata()
{
    this->calibration_data_names.push_back("tree");
    this->calibration_data_names.push_back("node_votes");
    this->calibration_data_names.push_back("decisions");

}


//modified
void calibration_model::read_calibration_models(const std::string& calibration_json_file,
                                                std::vector<std::string> variable_names)
{
    std::vector< set_of_calibrations_type > all_data;
    for (unsigned int v = 0; v < variable_names.size(); v++)
    {
       set_of_calibrations_type dummy_arr;
       all_data.push_back(dummy_arr);
    }
   try
   {
       boost::property_tree::ptree doc;
       std::string json_filename = calibration_json_file;
       boost::property_tree::read_json(json_filename, doc);

       BOOST_FOREACH(boost::property_tree::ptree::value_type &each_tree, doc)
       {
           for (unsigned int vn=0; vn < variable_names.size(); vn++)
           {

               std::map<int, std::vector<double> > node_votes;
               BOOST_FOREACH(boost::property_tree::ptree::value_type &v, each_tree.second.get_child(variable_names[vn]))
               {
                   std::vector<double> prob_tuple (2,0);
                   int ind = 0;
                   BOOST_FOREACH(boost::property_tree::ptree::value_type &i, v.second)
                   {
                       double p = i.second.get_value<double>();
                       prob_tuple[ind++] = p;
                   }
                   node_votes[atoi(v.first.c_str())] = prob_tuple;
               }
               all_data[vn].push_back(node_votes);
           }
       }
    int ind = 0;
    this->all_trees = all_data[ind++];
    this->all_node_votes = all_data[ind++];
    this->all_decisions = all_data[ind++];
   }

   catch (std::exception const& e)
   {
       std::cerr << e.what() << std::endl;
   }
}


double calibration_model::get_single_dectree_proba(const feature_type& features, int tree_index)
{
    int node = 0;
    while (this->all_trees[tree_index][node][0] != -1)
    {
        if (features.find((int)this->all_decisions[tree_index][node][0])->second <= this->all_decisions[tree_index][node][1])
        {
            node = (int)this->all_trees[tree_index][node][0];
        }
        else
        {
            node = (int)this->all_trees[tree_index][node][1];
        }
    }
    //node captures the index in the tree.

    double total = this->all_node_votes[tree_index][node][0] + this->all_node_votes[tree_index][node][1];
    double proba1 = this->all_node_votes[tree_index][node][1] / total;

    return proba1;

}

double calibration_model::get_randomforest_proba(const feature_type& features)
{
    //get the probability for every tree and average them out.
    double final_proba = 0;
    for (int t = 0; t < this->n_trees; t++)
    {
       final_proba += this->get_single_dectree_proba(features, t);
    }
    return final_proba/this->n_trees; // returns calibration score
}

double scoring_models::score_instance(const feature_type& features)
{
    return this->randomforest_model.get_randomforest_proba(features);
}

error_model& scoring_models::get_indel_model(const std::string& pattern){
    if (pattern=="f"){}
    return this->indel_models[this->current_indel_model].model;
}



void scoring_models::load_indel_models(boost::property_tree::ptree pt,const std::string& model_name){
    std::string s = imodels + "." + model_name;
//    log_os << s << std::endl;
    indel_model temp_model;
    unsigned i=0;
    BOOST_FOREACH(boost::property_tree::ptree::value_type &v, pt.get_child(s))
    {
        temp_model.add_prop(i,atof(v.second.data().c_str()),atof(v.second.data().c_str()));
        i++;
    }
    this->indel_models[model_name] = temp_model;
    this->indel_init = true;
}



void scoring_models::load_calibration_models(boost::property_tree::ptree pt,const std::string& model_name)
{
    this->randomforest_model.populate_storage_metadata();

    this->randomforest_model.read_calibration_models(calibration_json_file,
                                                     this->randomforest_model.calibration_data_names);

    this->calibration_init = true;
}

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
     this->load_calibration_models(calibration_json_file);

}

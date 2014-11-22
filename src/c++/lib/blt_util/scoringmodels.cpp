/*
 * scoringmodels.cpp
 *
 *  Created on: Aug 20, 2014
 *      Author: mkallberg
 */
#include "scoringmodels.hh"

//#define DEBUG_SCORINGMODELS

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
void calibration_model::load(boost::property_tree::ptree pt)
{
    std::vector< set_of_calibrations_type > all_data;
    for (unsigned int v = 0; v < this->calibration_data_names.size(); v++)
    {
       set_of_calibrations_type dummy_arr;
       all_data.push_back(dummy_arr);
    }
//   try
//   {
       int t_count = 0;

       BOOST_FOREACH(boost::property_tree::ptree::value_type &each_tree, pt)
       {
           t_count ++;
//           log_os << "Tree count: " << t_count << "\n";

           for (unsigned int vn=0; vn < this->calibration_data_names.size(); vn++)
           {

               std::map<int, std::vector<double> > node_votes;
               BOOST_FOREACH(boost::property_tree::ptree::value_type &v, each_tree.second.get_child(this->calibration_data_names[vn]))
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
//   }
//
//   catch (std::exception const& e)
//   {
//       std::cerr << e.what() << std::endl;
//   }
}


double calibration_model::get_single_dectree_proba(const feature_type& features, int tree_index)
{
    int node = 0;

    //traverse a single tree
    while (this->all_trees[tree_index][node][0] != -1)  // test condition signifies we've reached a leaf node
    {
//        log_os << "Looking for feature number " << STRELKA_VQSR_FEATURES::get_feature_label((int)this->all_decisions[tree_index][node][0]) << std::endl;
        if (features.find((int)this->all_decisions[tree_index][node][0])->second <= this->all_decisions[tree_index][node][1])
        {
//            log_os << "Looking for feature " << STRELKA_VQSR_FEATURES::get_feature_label((int)this->all_decisions[tree_index][node][0]) << std::endl;
            node = (int)this->all_trees[tree_index][node][0];
        }
        else
        {
            node = (int)this->all_trees[tree_index][node][1];
        }
    }

    // normalize the vote of the lead split
    double total = this->all_node_votes[tree_index][node][0] + this->all_node_votes[tree_index][node][1];
    double proba1 = 1.0*this->all_node_votes[tree_index][node][1] / total;

    return proba1;

}

double calibration_model::get_randomforest_proba(const feature_type& features)
{
    //get the probability for every tree and average them out.
    double final_proba = 0;
    for (int t = 0; t < this->n_trees; t++)
    {
//       log_os << "Applying tree " << t << std::endl;
       final_proba += this->get_single_dectree_proba(features, t);
    }
//    log_os << "Final prop " << 1-final_proba/this->n_trees << std::endl;
    return (1.0-final_proba/this->n_trees); // returns calibration score
}

int scoring_models::score_instance(const feature_type& features)
{
    double score = this->randomforest_model.get_randomforest_proba(features);
    int Q = error_prob_to_qphred(score);
    return Q;
}

const error_model& scoring_models::get_indel_model(const std::string& pattern){
    if (pattern=="f"){}
    return this->indel_models[this->current_indel_model].model;
}


void scoring_models::load_indel_model(boost::property_tree::ptree pt,const std::string& model_name){
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



void scoring_models::load_calibration_model(boost::property_tree::ptree pt,const std::string& model_name,const std::string& model_type)
{
    this->randomforest_model.populate_storage_metadata();
    if (model_name!="" && model_type!=""){
//        log_os << "Loading cali model: " << model_name << std::endl;
    }

    //TODO add case here for indels, currently only loading snps. Need to add another nesting level here
    this->randomforest_model.load(pt.get_child(model_name));
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
     BOOST_FOREACH(boost::property_tree::ptree::value_type &v, pt.get_child(imodels)){
//         log_os << "Reading indel model " << v.first <<  std::endl;
         this->load_indel_model(pt,v.first);
     }

     //load calibration models
     BOOST_FOREACH(boost::property_tree::ptree::value_type &v, pt.get_child(cmodels)){
         this->load_calibration_model(pt.get_child(cmodels),v.first);
     }

}

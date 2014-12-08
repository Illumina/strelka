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

#include "boost/property_tree/json_parser.hpp"

#include <cassert>
#include <cstdlib>     /* atof */

#include <iostream>
#include <sstream>

using boost::property_tree::ptree;

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

void calibration_model::populate_storage_metadata()
{
    this->calibration_data_names.push_back("tree");
    this->calibration_data_names.push_back("node_votes");
    this->calibration_data_names.push_back("decisions");

}


//modified
void calibration_model::load(const ptree& pt)
{
    const unsigned nameSize(calibration_data_names.size());
    std::vector< set_of_calibrations_type > all_data(nameSize);

    //   try
//   {
    int t_count = 0;

    for (const ptree::value_type& each_tree : pt)
    {
        t_count ++;
//           log_os << "Tree count: " << t_count << "\n";

        for (unsigned int vn=0; vn < nameSize; vn++)
        {

            std::map<int, std::vector<double> > node_votes;
            for (const ptree::value_type& v : each_tree.second.get_child(this->calibration_data_names[vn]))
            {
                std::vector<double> prob_tuple (2,0);
                int ind = 0;
                for (const ptree::value_type& i : v.second)
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


double calibration_model::get_single_dectree_proba(const feature_type& features, int tree_index) const
{
    const calibration_type& tree(all_trees[tree_index]);
    const calibration_type& decision(all_decisions[tree_index]);
    int node = 0;

    //traverse a single tree
    while (tree.at(node)[0] != -1)  // test condition signifies we've reached a leaf node
    {
//        log_os << "Looking for feature number " << STRELKA_VQSR_FEATURES::get_feature_label((int)this->all_decisions[tree_index][node][0]) << std::endl;
        if (features.at((int)decision.at(node)[0]) <= decision.at(node)[1])
        {
//            log_os << "Looking for feature " << STRELKA_VQSR_FEATURES::get_feature_label((int)this->all_decisions[tree_index][node][0]) << std::endl;
            node = (int)tree.at(node)[0];
        }
        else
        {
            node = (int)tree.at(node)[1];
        }
    }

    // normalize the vote of the lead split
    const std::vector<double>& votes(all_node_votes[tree_index].at(node));
    const double total = votes[0] + votes[1];
    return (votes[1] / total);
}



double calibration_model::get_randomforest_proba(const feature_type& features) const
{
    double retval(0);
    try
    {
        //get the probability for every tree and average them out.
        double final_proba = 0;
        for (int t = 0; t < this->n_trees; t++)
        {
            final_proba += this->get_single_dectree_proba(features, t);
        }
        retval = (1.0-final_proba/this->n_trees);
    }
    catch (...)
    {
        log_os << "Except caught in random forest while scoring feature:\n";
        for (const auto val : features)
        {
            log_os << "K:V " << val.first << " : " << val.second << "\n";
        }
        throw;
    }
    return retval;
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
    const double score = this->randomforest_model.get_randomforest_proba(features);
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
    for (const ptree::value_type &v : pt.get_child(s))
    {
        temp_model.add_prop(i,atof(v.second.data().c_str()),atof(v.second.data().c_str()));
        i++;
    }
    this->indel_models[model_name] = temp_model;
    this->indel_init = true;
}



void scoring_models::load_calibration_model(const ptree& pt,const std::string& model_name,const std::string& model_type)
{
    this->randomforest_model.populate_storage_metadata();
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
    for (const ptree::value_type &v : pt.get_child(imodels))
    {
//         log_os << "Reading indel model " << v.first <<  std::endl;
        this->load_indel_model(pt,v.first);
    }

    //load calibration models
    for (const ptree::value_type &v : pt.get_child(cmodels))
    {
        this->load_calibration_model(pt.get_child(cmodels),v.first);
    }
}

// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2016 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//
/*
 * scoringmodels.cpp
 *
 *  Created on: Aug 20, 2014
 *      Author: Morten Kallberg
 */

#ifdef _MSC_VER
#define NOMINMAX
#endif

#include "scoringmodels.hh"

#include "blt_util/log.hh"
#include "blt_util/qscore.hh"
#include "blt_util/blt_exception.hh"

#include <cassert>
#include <fstream>
#include <algorithm>
#include <iostream>

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
scoring_models& scoring_models::Instance()
{
    if (!m_pInstance)    // Only allow one instance of class to be generated.
    {
        m_pInstance = new scoring_models;
    }
    return *m_pInstance;
}

void scoring_models::init_default_models()
{
    // load previously hard-coded logistic regression model
    IndelErrorModel newModel = generate_new_indel_error_model();
    this->indel_models[newModel.getName()] = newModel;

    // load previously hard-coded polynomial model
    IndelErrorModel oldModel = generate_old_indel_error_model();
    this->indel_models[oldModel.getName()] = oldModel;
}



const IndelErrorModel& scoring_models::get_indel_model() const
{
    assert(!this->current_indel_model.empty());
    return this->indel_models.at(this->current_indel_model);
}


// Generated the header VCF header string specifying which json file and what models were used
void
scoring_models::
writeVcfHeader(
    std::ostream& os) const
{
    if (! indelErrorModelFilename.empty())
    {
        os << "##IndelErrorModelFilename=" << this->indelErrorModelFilename << "\n";
    }

    if (! current_indel_model.empty())
    {
        os << "##IndelModel=" << this->current_indel_model << "\n";
    }
}


void scoring_models::set_indel_model(const std::string& model_name)
{
    // if indel error model is already set, it cannot be reset
    assert(this->current_indel_model.empty() && "Cannot reset indel error model once set");

    // check to see that indel error model key is in model map
    if (indel_models.find(model_name) == indel_models.end())
    {
        std::string model_names;
        possible_indel_models(model_names);
        std::string errmsg = "Cannot find " + model_name + " in indel model map\n";
        errmsg += "Possible indel models are: " + model_names;
        throw blt_exception(errmsg.c_str());
    }
    current_indel_model = model_name;
}

void scoring_models::possible_indel_models(std::string& str) const
{

    if (this->indel_models.empty())
    {
        str = "No indel models loaded";
        return;
    }

    std::map<std::string,IndelErrorModel>::const_iterator it = indel_models.begin();
    str = it->first;
    for (; it != indel_models.end(); ++it)
        str += "," + it->first;
}

void scoring_models::load_indel_model(const Json::Value& data)
{
    IndelErrorModel tempModel;
    tempModel.Deserialize(data);
    this->indel_models[tempModel.getName()] = tempModel;
}



void
scoring_models::
load_single_model(
    const Json::Value& data,
    const MODEL_TYPE::index_t model_type)
{
    using namespace MODEL_TYPE;
    switch (model_type)
    {
    case INDELMODEL:
        this->load_indel_model(data);
        break;
    default:
        assert(false && "Unknown model-type in model json.");
    }
}



static
void
loadFileToJson(
    const std::string& filename,
    Json::Value& root)
{
    std::ifstream file( filename , std::ifstream::binary);
    file >> root;
}



void
scoring_models::
load_models(
    const std::string& model_file,
    const MODEL_TYPE::index_t model_type)
{
    Json::Value root;
    loadFileToJson(model_file,root);

    Json::Value models = root[get_label(model_type)];
    if (models.isNull()) return;

    for (auto& modelValue : models)
    {
        load_single_model(modelValue, model_type);
    }
}



void
scoring_models::
load_indel_error_models(
    const std::string& model_file)
{
    indelErrorModelFilename = model_file;
    load_models(model_file,MODEL_TYPE::INDELMODEL);
}

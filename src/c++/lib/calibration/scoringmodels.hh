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
 * scoringmodels.hh
 *
 *  Created on: Aug 20, 2014
 *      Author: Morten Kallberg
 */

#pragma once

#include "calibration/Indelmodel.hh"

#include "json/json.h"

#include <iosfwd>
#include <map>
#include <string>


namespace MODEL_TYPE
{
enum index_t
{
    INDELMODEL,
    SIZE
};

inline
const std::string
get_label(const index_t i)
{
    switch (i)
    {
    case INDELMODEL:
        return "IndelModels";
    default:
        assert(false && "Unknown node type in scoringmodels");
        return nullptr;
    }
}
}

struct scoring_models
{
    static scoring_models& Instance();

    void load_indel_error_models(const std::string& model_file);

    void set_indel_model(const std::string& model_name);

    void get_indel_error(
        const starling_base_options& client_opt,
        const starling_indel_report_info& iri,
        double& indel_error_prob,
        double& ref_error_prob) const
    {
        get_indel_model().calc_prop(client_opt,iri,indel_error_prob,ref_error_prob);
    }

    const IndelErrorModel& get_indel_model() const;

    bool isIndelInit() const
    {
        return (! indel_models.empty());
    }

    void
    writeVcfHeader(
        std::ostream& os) const;
private:
    scoring_models()
    {
        this->init_default_models();
    } // Private so that it cannot be called
    scoring_models(scoring_models const&) {}            // copy constructor is private
    scoring_models& operator=(scoring_models const&);  // assignment operator is private

    void init_default_models();

    void
    load_models(
        const std::string& model_file,
        const MODEL_TYPE::index_t model_type);

    void
    load_single_model(
        const Json::Value& data,
        const MODEL_TYPE::index_t model_type);

    void load_indel_model(const Json::Value& data);

    const IndelErrorModel& get_indel_model(const std::string& pattern) const;

    void possible_indel_models(std::string& str) const;

    static scoring_models* m_pInstance;

    std::string indelErrorModelFilename;

    typedef std::map<std::string,IndelErrorModel> indel_modelmap;
    indel_modelmap indel_models;
    std::string current_indel_model;
};

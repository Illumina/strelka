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
 *      Author: Morten Kallberg
 */

#pragma once

#include "calibration/RandomForestModel.hh"
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
    CALMODEL,
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
    case CALMODEL:
        return "CalibrationModels";
    default:
        assert(false && "Unknown node type in scoringmodels");
        return nullptr;
    }
}
}

struct scoring_models
{
    static scoring_models& Instance();

    void load_variant_scoring_models(const std::string& model_file);
    void load_indel_error_models(const std::string& model_file);

    void set_indel_model(const std::string& model_name);

    double score_variant(
        const feature_type& features,
        const VARIATION_NODE_TYPE::index_t vtype) const;

    void get_indel_error(
        const starling_base_options& client_opt,
        const starling_indel_report_info& iri,
        double& indel_error_prob,
        double& ref_error_prob) const
    {
        get_indel_model().calc_prop(client_opt,iri,indel_error_prob,ref_error_prob);
    }

    const IndelErrorModel& get_indel_model() const;

    void writeVcfHeader(std::ostream& os) const;

    bool isVariantScoringInit() const
    {
        return randomforest_model.isInit();
    }
    bool isIndelInit() const
    {
        return (! indel_models.empty());
    }

    // See STARKA-257 github comment for more detail on this fit
    double recal_somatic_snv_score(double& score) const {
    	return 2.57*score+0.94;
    }

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
    void load_calibration_model(const Json::Value& data);

    const IndelErrorModel& get_indel_model(const std::string& pattern) const;

    void possible_indel_models(std::string& str) const;

    static scoring_models* m_pInstance;

    std::string variantScoringModelFilename;
    std::string indelErrorModelFilename;

    typedef std::map<std::string,IndelErrorModel> indel_modelmap;
    indel_modelmap indel_models;
    std::string current_indel_model;

    RandomForestModel randomforest_model,randomforest_model_indel;
};

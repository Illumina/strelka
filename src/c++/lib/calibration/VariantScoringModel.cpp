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

#include "VariantScoringModel.hh"
#include "blt_util/blt_exception.hh"

#include <fstream>
#include <iostream>
#include <sstream>



static
void
modelParseError(
    const std::string& model_file,
    const std::string& key)
{
    std::ostringstream oss;
    oss << "Can't find node '" << key << "' in scoring model file: '" << model_file << "'";
    throw blt_exception(oss.str().c_str());
}



VariantScoringModel::
VariantScoringModel(
    const std::string& model_file,
    const SCORING_CALL_TYPE::index_t ctype,
    const SCORING_VARIANT_TYPE::index_t vtype)
{
    Json::Value root;
    {
        std::ifstream file(model_file , std::ifstream::binary);
        file >> root;
    }

    static const std::string model_type("CalibrationModels");
    Json::Value models = root[model_type];
    if (models.isNull()) modelParseError(model_file,model_type);

    const std::string call_type(SCORING_CALL_TYPE::get_label(ctype));
    Json::Value callmodels = models[call_type];
    if (callmodels.isNull()) modelParseError(model_file,call_type);

    const std::string var_type(SCORING_VARIANT_TYPE::get_label(vtype));
    Json::Value varmodel = callmodels[var_type];
    if (varmodel.isNull()) modelParseError(model_file,var_type);

    _meta.Deserialize(varmodel);
    _model.Deserialize(varmodel);
}

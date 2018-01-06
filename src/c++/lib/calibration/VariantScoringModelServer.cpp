//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2018 Illumina, Inc.
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
 *      Author: Morten Kallberg
 */

#include "VariantScoringModelServer.hh"

#include "RandomForestModel.hh"

#include "blt_util/log.hh"
#include "common/Exceptions.hh"

#include "rapidjson/filereadstream.h"

#include <cstdio>

#include <fstream>
#include <iostream>
#include <sstream>



static
void
modelParseError(
    const std::string& modelFile,
    const std::string& key)
{
    std::ostringstream oss;
    oss << "Can't find node '" << key << "' in json scoring model file: '" << modelFile << "'";
    BOOST_THROW_EXCEPTION(illumina::common::GeneralException(oss.str()));
}



VariantScoringModelServer::
VariantScoringModelServer(
    const VariantScoringModelMetadata::featureMap_t& featureMap,
    const std::string& modelFile,
    const SCORING_CALL_TYPE::index_t callType,
    const SCORING_VARIANT_TYPE::index_t variantType)
{
    rapidjson::Document document;
    {
        FILE* tmpFilePtr = fopen(modelFile.c_str(), "rb");
        char readBuffer[65536];
        rapidjson::FileReadStream inputFileStream(tmpFilePtr, readBuffer, sizeof(readBuffer));
        if (document.ParseStream(inputFileStream).HasParseError())
        {
            std::ostringstream oss;
            oss << "ERROR: Failed to parse json scoring model file: '" << modelFile << "'";
            BOOST_THROW_EXCEPTION(illumina::common::GeneralException(oss.str()));
        }
        fclose(tmpFilePtr);
    }

    if (! document.IsObject())
    {
        std::ostringstream oss;
        oss << "Unexpected root data type in json scoring model file: '" << modelFile << "'";
        BOOST_THROW_EXCEPTION(illumina::common::GeneralException(oss.str()));
    }

    auto getNodeMember = [&](const rapidjson::Value& node, const char* label) -> const rapidjson::Value&
    {
        const rapidjson::Value::ConstMemberIterator iter(node.FindMember(label));
        if (iter == node.MemberEnd()) modelParseError(modelFile, label);
        return iter->value;
    };

    static const std::string modelTypeLabel("CalibrationModels");
    const rapidjson::Value& models(getNodeMember(document, modelTypeLabel.c_str()));

    const rapidjson::Value& callModels(getNodeMember(models, SCORING_CALL_TYPE::get_label(callType)));

    const rapidjson::Value& varModels(getNodeMember(callModels, SCORING_VARIANT_TYPE::get_label(variantType)));

    try
    {
        _meta.Deserialize(featureMap, varModels);

        if (_meta.modelType == "RandomForest")
        {
            std::unique_ptr<RandomForestModel> rfModel(new RandomForestModel());
            rfModel->Deserialize(featureMap.size(), varModels);
            _model = std::move(rfModel);
        }
        else
        {
            assert(false && "Unrecognized model type");
        }
    }
    catch (...)
    {
        log_os << "ERROR: Exception caught while attempting to parse scoring model file '" << modelFile << "'\n";
        throw;
    }
}

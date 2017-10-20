//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2017 Illumina, Inc.
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

#include <iomanip>
#include <iostream>
#include <fstream>

#include "rapidjson/filereadstream.h"
#include "common/Exceptions.hh"
#include "blt_util/log.hh"
#include "IndelErrorModelJson.hh"
#include "ThetaJson.hh"
#include "common/RapidJsonHelper.hh"

IndelErrorModelJson::
IndelErrorModelJson(
        const std::string& sampleName,
        const IndelErrorModelBinomialMixture& model,
        const bool isStatic)
    : _sampleName(sampleName)
    , _model(model)
    , _isStatic(isStatic)
{}

IndelErrorModelJson
IndelErrorModelJson::deserialize(
        const rapidjson::Value& root)
{
    using namespace illumina::common;

    static const char* sampleNameLabel = "sampleName";
    const rapidjson::Value& sampleNameValue(RapidJsonHelper::getNodeMember(root, sampleNameLabel));

    const std::string sampleName(sampleNameValue.GetString());

    static const char* motifLabel = "motif";
    const rapidjson::Value& motifArray(RapidJsonHelper::getNodeMember(root, motifLabel));

    IndelErrorModelBinomialMixture model(IndelErrorModelBinomialMixture::deserialize(motifArray));

    static const char* isStaticLabel = "isStatic";
    const rapidjson::Value& isStaticValue(RapidJsonHelper::getNodeMember(root, isStaticLabel));

    const bool isStatic(isStaticValue.GetBool());

    return IndelErrorModelJson(sampleName, model, isStatic);
}

IndelErrorModelsJson
IndelErrorModelsJson::deserialize(
        const rapidjson::Value& root)
{
    using namespace illumina::common;

    IndelErrorModelsJson indelErrorModelsJson;
    static const char* sampleLabel = "sample";
    const rapidjson::Value& sampleArray(RapidJsonHelper::getNodeMember(root, sampleLabel));

    // one json file could potentially have multiple samples
    for (const auto& sampleValue : sampleArray.GetArray())
    {
        indelErrorModelsJson.addModel(IndelErrorModelJson::deserialize(sampleValue));
    }

    return indelErrorModelsJson;
}

void
IndelErrorModelParser::importIndelErrorModelJsonFile(
    const std::string& modelFilename,
    std::map<std::string, IndelErrorRateSet>& modelMap)
{

    IndelErrorModelsJson indelErrorModelsJson;
    importIndelErrorModelJsonFile(modelFilename, indelErrorModelsJson);

    const auto indelErrorModels(indelErrorModelsJson.getIndelErrorModels());

    for (const auto& sampleIndelErrorModel : indelErrorModels)
    {
        IndelErrorRateSet& modelSample(modelMap[sampleIndelErrorModel.getSampleName()]);
        const auto binomialMixtureModelMotifs(sampleIndelErrorModel.getBinomialMixtureModel().getMotifs());

        for (const auto& motif : binomialMixtureModelMotifs)
        {
            modelSample.addRate(motif.getRepeatPatternSize(), motif.getRepeatCount(), motif.getIndelRate(), motif.getIndelRate(), motif.getNoisyLocusRate());
        }
    }
}

void
IndelErrorModelParser::importIndelErrorModelJsonFile(
    const std::string& modelFilename,
    IndelErrorModelsJson& indelErrorModelsJson)
{
    rapidjson::Document document;
    {
        FILE* tmpFilePtr = fopen(modelFilename.c_str(), "rb");
        char readBuffer[65536];
        rapidjson::FileReadStream inputFileStream(tmpFilePtr, readBuffer, sizeof(readBuffer));
        if (document.ParseStream(inputFileStream).HasParseError())
        {
            std::ostringstream oss;
            oss << "ERROR: Failed to parse json indel parameter file: '" << modelFilename << "'";
            BOOST_THROW_EXCEPTION(illumina::common::LogicException(oss.str()));
        }
        fclose(tmpFilePtr);
    }

    if (!document.IsObject())
    {
        std::ostringstream oss;
        oss << "ERROR: Unexpected root data type in json indel parameter file: '" << modelFilename << "'";
        BOOST_THROW_EXCEPTION(illumina::common::LogicException(oss.str()));
    }

    indelErrorModelsJson = IndelErrorModelsJson::deserialize(document);
}

std::map<std::string, IndelErrorRateSet>
IndelErrorModelParser::
generateIndelErrorRateSetMap(
        const std::vector<std::string>& modelFilenames)
{
    std::map<std::string, IndelErrorRateSet> modelMap;
    for (const auto& modelFilename : modelFilenames)
    {
        importIndelErrorModelJsonFile(modelFilename, modelMap);
    }
    return modelMap;
}

void
IndelErrorModelParser::
importThetaJsonFile(
    const std::string& thetaFilename,
    std::map<unsigned, std::vector<double>>& thetasMap)
{
    rapidjson::Document document;
    {
        FILE* tmpFilePtr = fopen(thetaFilename.c_str(), "rb");
        char readBuffer[65536];
        rapidjson::FileReadStream inputFileStream(tmpFilePtr, readBuffer, sizeof(readBuffer));
        if (document.ParseStream(inputFileStream).HasParseError())
        {
            std::ostringstream oss;
            oss << "ERROR: Failed to parse json theta file: '" << thetaFilename << "'";
            BOOST_THROW_EXCEPTION(illumina::common::LogicException(oss.str()));
        }
        fclose(tmpFilePtr);
    }

    if (!document.IsObject())
    {
        std::ostringstream oss;
        oss << "ERROR: Unexpected root data type in json theta file: '" << thetaFilename << "'";
        BOOST_THROW_EXCEPTION(illumina::common::LogicException(oss.str()));
    }

    thetasMap = ThetasJson::deserialize(document).getThetasMap();
}


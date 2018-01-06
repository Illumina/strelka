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
            oss << "Failed to parse json indel parameter file: '" << modelFilename << "'";
            BOOST_THROW_EXCEPTION(illumina::common::GeneralException(oss.str()));
        }
        fclose(tmpFilePtr);
    }

    if (!document.IsObject())
    {
        std::ostringstream oss;
        oss << "ERROR: Unexpected root data type in json indel parameter file: '" << modelFilename << "'";
        BOOST_THROW_EXCEPTION(illumina::common::GeneralException(oss.str()));
    }

    try
    {
        indelErrorModelsJson = IndelErrorModelsJson::deserialize(document);
    }
    catch (...)
    {
        log_os << "Exception caught while deserializing json file '" << modelFilename << "'\n";
        throw;
    }
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
            oss << "Failed to parse json theta file: '" << thetaFilename << "'";
            BOOST_THROW_EXCEPTION(illumina::common::GeneralException(oss.str()));
        }
        fclose(tmpFilePtr);
    }

    if (!document.IsObject())
    {
        std::ostringstream oss;
        oss << "Unexpected root data type in json theta file: '" << thetaFilename << "'";
        BOOST_THROW_EXCEPTION(illumina::common::GeneralException(oss.str()));
    }

    try
    {
        thetasMap = ThetasJson::deserialize(document).getThetasMap();
    }
    catch (...)
    {
        log_os << "Exception caught while deserializing json file '" << thetaFilename << "'\n";
        throw;
    }
}


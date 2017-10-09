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

#include "rapidjson/istreamwrapper.h"
#include "common/Exceptions.hh"
#include "blt_util/log.hh"
#include "IndelErrorModelJson.hh"

IndelErrorModelJson::
IndelErrorModelJson(const std::string& sampleName, const IndelErrorModelBinomialMixture& model, const bool isStatic)
    : _sampleName(sampleName)
    , _model(model)
    , _isStatic(isStatic)
{}


std::map<std::string, IndelErrorRateSet>
IndelErrorModelJson::
generateIndelErrorRateSetMap(const std::vector<std::string>& modelFilenames)
{
    std::map<std::string, IndelErrorRateSet> modelMap;
    for (const auto& modelFilename : modelFilenames)
    {
        loadIndelErrorRateSet(modelFilename, modelMap);
    }
    return modelMap;
}

void
IndelErrorModelJson::
loadIndelErrorRateSet(
    const std::string& modelFilename,
    std::map<std::string, IndelErrorRateSet>& modelMap)
{
    rapidjson::Document document;
    std::ifstream ifs(modelFilename);

    if (!ifs.is_open())
    {
        using namespace illumina::common;
        std::ostringstream oss;
        oss << "ERROR: Cannot open indel error model file '" << modelFilename << "'\n";
        BOOST_THROW_EXCEPTION(LogicException(oss.str()));
    }
    rapidjson::IStreamWrapper isw(ifs);
    document.ParseStream(isw);

    const rapidjson::Value& samples(document["sample"]);
    if (samples.Empty() || !samples.IsArray())
    {
        using namespace illumina::common;
        std::ostringstream oss;
        oss << "ERROR: no samples in indel error model file '" << modelFilename << "'\n";
        BOOST_THROW_EXCEPTION(LogicException(oss.str()));
    }

    // one json file could potentially have multiple samples
    for (rapidjson::Value::ConstValueIterator sample = samples.Begin(); sample != samples.End(); ++sample)
    {
        const std::string sampleName((*sample)["sampleName"].GetString());
        modelMap[sampleName] = IndelErrorRateSet();
        const rapidjson::Value& motifs((*sample)["motif"]);
        if (motifs.Empty()|| !motifs.IsArray())
        {
            using namespace illumina::common;
            std::ostringstream oss;
            oss << "ERROR: no params for sample '" << sampleName << "' in indel error model file '" << modelFilename << "'\n";
            BOOST_THROW_EXCEPTION(LogicException(oss.str()));
        }

        for (rapidjson::Value::ConstValueIterator motif = motifs.Begin(); motif != motifs.End(); ++motif)
        {
            const double indelRate = (*motif)["indelRate"].GetDouble();
            const double noisyLocusRate = (*motif)["noisyLocusRate"].GetDouble();
            const unsigned repeatCount = (*motif)["repeatCount"].GetUint();
            const unsigned repeatPatternSize = (*motif)["repeatPatternSize"].GetUint();
            modelMap[sampleName].addRate(repeatPatternSize, repeatCount, indelRate, indelRate, noisyLocusRate);
        }
    }
}


std::map<unsigned, std::vector<double> >
IndelErrorModelJson::
deserializeTheta(
    const std::string& filename)
{
    std::map<unsigned, std::vector<double>> thetasMap;

    rapidjson::Document document;

    std::ifstream ifs(filename, std::ifstream::binary);

    if (!ifs.is_open())
    {
        using namespace illumina::common;
        std::ostringstream oss;
        oss << "ERROR: Cannot open theta file '" << filename << "'\n";
        BOOST_THROW_EXCEPTION(LogicException(oss.str()));
    }

    rapidjson::IStreamWrapper isw(ifs);
    document.ParseStream(isw);

    const rapidjson::Value& thetas(document["thetas"]);
    if (thetas.Empty() || !thetas.IsArray())
    {
        using namespace illumina::common;
        std::ostringstream oss;
        oss << "ERROR: no thetas in theta file '" << filename << "'\n";
        BOOST_THROW_EXCEPTION(LogicException(oss.str()));
    }

    for (rapidjson::Value::ConstValueIterator theta = thetas.Begin(); theta != thetas.End(); ++theta)
    {
        std::vector<double> thetaVector;
        unsigned repeatPatternSize = (*theta)["repeatPatternSize"].GetUint();
        const rapidjson::Value& thetaValues((*theta)["theta"]);

        if (thetaValues.Empty() || !thetaValues.IsArray())
        {
            using namespace illumina::common;
            std::ostringstream oss;
            oss << "ERROR: no theta values in theta file '" << filename << "'\n";
            BOOST_THROW_EXCEPTION(LogicException(oss.str()));
        }

        for (rapidjson::Value::ConstValueIterator thetaValue = thetaValues.Begin(); thetaValue != thetaValues.End(); ++thetaValue)
        {
            thetaVector.push_back(thetaValue->GetDouble());
        }
        thetasMap[repeatPatternSize] = thetaVector;
    }
    return thetasMap;
}

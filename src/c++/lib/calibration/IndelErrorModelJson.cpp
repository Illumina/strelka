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

#include "common/Exceptions.hh"
#include "blt_util/log.hh"
#include "IndelErrorModelJson.hh"

IndelErrorModelJson::
IndelErrorModelJson(const std::string& sampleName)
    : _sampleName(sampleName)
{}


// move these to a more appropriate place later
Json::Value
IndelErrorModelJson::
generateMotifsNode() const
{
    Json::Value motifs;
    for (const auto& motifIt : model.motifs)
    {
        Json::Value motif;
        motif["repeatPatternSize"] = motifIt.repeatPatternSize;
        motif["repeatCount"] = motifIt.repeatCount;
        motif["indelRate"] = motifIt.indelRate;
        motif["noisyLocusRate"] = motifIt.noisyLocusRate;
        motifs.append(motif);
    }
    return motifs;
}

void
IndelErrorModelJson::
serializeIndelErrorModel(
    const std::string& sampleName,
    const Json::Value& motifsNode,
    const bool isStatic,
    const std::string& filename)
{
    Json::StyledWriter writer;
    Json::Value jsonRoot;
    Json::Value samples;
    Json::Value sample;
    sample["sampleName"] = sampleName;
    sample["motif"] = motifsNode;
    sample["isStatic"] = isStatic;
    samples.append(sample);
    jsonRoot["sample"] = samples;

    const std::string str = writer.write(jsonRoot);
    std::ofstream out(filename);
    out << str << "\n\n";
}

void
IndelErrorModelJson::
addMotif(unsigned repeatPatternSize,
         unsigned repeatCount,
         double indelRate,
         double noisyLocusRate)
{
    IndelMotifBinomialMixture motif;
    motif.repeatPatternSize = repeatPatternSize;
    motif.repeatCount = repeatCount;
    motif.indelRate = indelRate;
    motif.noisyLocusRate = noisyLocusRate;
    model.motifs.push_back(motif);
}

std::map<std::string, IndelErrorRateSet>
IndelErrorModelJson::
deserializeIndelErrorModels(const std::vector<std::string>& modelFilenames)
{
    std::map<std::string, IndelErrorRateSet> modelMap;
    for (const auto& modelFilename : modelFilenames)
    {
        std::string jsonString;
        Json::Value root;
        {
            std::ifstream ifs(modelFilename, std::ifstream::binary);
            std::stringstream buffer;
            buffer << ifs.rdbuf();
            jsonString = buffer.str();
        }
        Json::Reader reader;
        if (!reader.parse(jsonString, root))
        {
            using namespace illumina::common;

            std::ostringstream oss;
            oss << "Failed to parse JSON " << modelFilename << " " << reader.getFormattedErrorMessages() << "'\n";
            BOOST_THROW_EXCEPTION(LogicException(oss.str()));
        }

        Json::Value samples = root["sample"];
        if (samples.isNull() || samples.empty())
        {
            using namespace illumina::common;
            std::ostringstream oss;
            oss << "ERROR: no samples in indel error model file '" << modelFilename << "'\n";
            BOOST_THROW_EXCEPTION(LogicException(oss.str()));
        }

        // one json file could potentially have multiple samples
        for (const auto& sample : samples)
        {
            std::string sampleName = sample["sampleName"].asString();
            modelMap[sampleName] = IndelErrorRateSet();
            Json::Value motifs = sample["motif"];
            if (motifs.isNull() || motifs.empty())
            {
                using namespace illumina::common;
                std::ostringstream oss;
                oss << "ERROR: no params for sample '" << sampleName << "' in indel error model file '" << modelFilename << "'\n";
                BOOST_THROW_EXCEPTION(LogicException(oss.str()));
            }

            for (const auto& motifValue : motifs)
            {
                const double indelRate = motifValue["indelRate"].asDouble();
                const double noisyLocusRate = motifValue["noisyLocusRate"].asDouble();
                const unsigned repeatCount = motifValue["repeatCount"].asInt();
                const unsigned repeatPatternSize = motifValue["repeatPatternSize"].asInt();
                modelMap[sampleName].addRate(repeatPatternSize, repeatCount, indelRate, indelRate, noisyLocusRate);
            }
        }

    }
    return modelMap;
}

std::map<unsigned, std::vector<double> >
IndelErrorModelJson::
deserializeTheta(
    const std::string& filename)
{
    std::string jsonString;
    Json::Value root;
    {
        std::ifstream ifs(filename, std::ifstream::binary);
        std::stringstream buffer;
        buffer << ifs.rdbuf();
        jsonString = buffer.str();
    }
    Json::Reader reader;
    reader.parse(jsonString, root);
    Json::Value thetasRoot = root["thetas"];
    if (thetasRoot.isNull() || thetasRoot.empty())
    {
        using namespace illumina::common;
        std::ostringstream oss;
        oss << "ERROR: no theta values in theta file '" << filename << "'\n";
        BOOST_THROW_EXCEPTION(LogicException(oss.str()));
    }

    std::map<unsigned, std::vector<double>> thetas;

    for (const auto& thetasByPatternSize : thetasRoot)
    {
        std::vector<double> theta;
        unsigned repeatPatternSize = thetasByPatternSize["repeatPatternSize"].asUInt();
        Json::Value thetaValues = thetasByPatternSize["theta"];
        for (const auto& thetaValue : thetaValues)
        {
            theta.push_back(thetaValue.asDouble());
        }
        thetas[repeatPatternSize] = theta;
    }
    return thetas;
}

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

#include "rapidjson/document.h"
#include "rapidjson/filereadstream.h"
#include "common/Exceptions.hh"
#include "blt_util/log.hh"
#include "IndelErrorModelJson.hh"



IndelErrorModelJson::
IndelErrorModelJson(const std::string& sampleName, const IndelErrorModelBinomialMixture& model, const bool isStatic)
    : _sampleName(sampleName)
    , _model(model)
    , _isStatic(isStatic)
{}



static
void
missingNodeError(
    const std::string& fileName,
    const char* key)
{
    std::ostringstream oss;
    oss << "ERROR: Can't find expected node '" << key << "' in json indel parameter file '" << fileName << "'.";
    BOOST_THROW_EXCEPTION(illumina::common::LogicException(oss.str()));
}



static
void
wrongValueTypeError(
    const std::string& fileName,
    const char* key,
    const char* keyType)
{
    std::ostringstream oss;
    oss << "ERROR: Node '" << key << "' does not have expected type '" << keyType
        << "' in json indel parameter file '" << fileName << "'.";
    BOOST_THROW_EXCEPTION(illumina::common::LogicException(oss.str()));
}



static
void
deserializeIndelErrorRateSet(
    const std::string& modelFilename,
    std::map<std::string, IndelErrorRateSet>& modelMap)
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

    auto getNodeMember = [&](
                             const rapidjson::Value& node,
                             const char* label) -> const rapidjson::Value&
    {
        const rapidjson::Value::ConstMemberIterator iter(node.FindMember(label));
        if (iter == node.MemberEnd()) missingNodeError(modelFilename, label);
        return iter->value;
    };

    static const char* sampleLabel = "sample";
    const rapidjson::Value& sampleArray(getNodeMember(document, sampleLabel));
    if (!sampleArray.IsArray()) wrongValueTypeError(modelFilename, sampleLabel, "array");

    if (sampleArray.Empty())
    {
        using namespace illumina::common;
        std::ostringstream oss;
        oss << "ERROR: No samples in json indel parameter file '" << modelFilename << "'\n";
        BOOST_THROW_EXCEPTION(LogicException(oss.str()));
    }

    // one json file could potentially have multiple samples
    for (const auto& sampleValue : sampleArray.GetArray())
    {
        static const char* sampleNameLabel = "sampleName";
        const rapidjson::Value& sampleNameValue(getNodeMember(sampleValue, sampleNameLabel));
        if (!sampleNameValue.IsString()) wrongValueTypeError(modelFilename, sampleNameLabel, "string");

        const std::string sampleName(sampleNameValue.GetString());

        static const char* motifLabel = "motif";
        const rapidjson::Value& motifArray(getNodeMember(sampleValue, motifLabel));
        if (!motifArray.IsArray()) wrongValueTypeError(modelFilename, motifLabel, "array");

        if (motifArray.Empty())
        {
            std::ostringstream oss;
            oss << "ERROR: No indel error rate parameters for sample '" << sampleName
                << "' in indel parameter file '" << modelFilename << "'\n";
            BOOST_THROW_EXCEPTION(illumina::common::LogicException(oss.str()));
        }

        IndelErrorRateSet& modelSample(modelMap[sampleName]);
        for (const auto& motifValue : motifArray.GetArray())
        {
            static const char* indelRateLabel = "indelRate";
            const rapidjson::Value& indelRateValue(getNodeMember(motifValue, indelRateLabel));
            if (!indelRateValue.IsNumber()) wrongValueTypeError(modelFilename, indelRateLabel, "number");
            const double indelRate(indelRateValue.GetDouble());

            static const char* noisyLocusRateLabel = "noisyLocusRate";
            const rapidjson::Value& noisyLocusRateValue(getNodeMember(motifValue, noisyLocusRateLabel));
            if (!noisyLocusRateValue.IsNumber()) wrongValueTypeError(modelFilename, noisyLocusRateLabel, "number");
            const double noisyLocusRate(noisyLocusRateValue.GetDouble());

            static const char* repeatCountLabel = "repeatCount";
            const rapidjson::Value& repeatCountValue(getNodeMember(motifValue, repeatCountLabel));
            if (!repeatCountValue.IsUint()) wrongValueTypeError(modelFilename, repeatCountLabel, "unsigned");
            const unsigned repeatCount(repeatCountValue.GetUint());

            static const char* repeatPatternSizeLabel = "repeatPatternSize";
            const rapidjson::Value& repeatPatternSizeValue(getNodeMember(motifValue, repeatPatternSizeLabel));
            if (!repeatPatternSizeValue.IsUint())
                wrongValueTypeError(modelFilename, repeatPatternSizeLabel, "unsigned");
            const unsigned repeatPatternSize(repeatPatternSizeValue.GetUint());

            modelSample.addRate(repeatPatternSize, repeatCount, indelRate, indelRate, noisyLocusRate);
        }
    }
}



std::map<std::string, IndelErrorRateSet>
IndelErrorModelJson::
generateIndelErrorRateSetMap(const std::vector<std::string>& modelFilenames)
{
    std::map<std::string, IndelErrorRateSet> modelMap;
    for (const auto& modelFilename : modelFilenames)
    {
        deserializeIndelErrorRateSet(modelFilename, modelMap);
    }
    return modelMap;
}



std::map<unsigned, std::vector<double> >
IndelErrorModelJson::
deserializeTheta(
    const std::string& thetaFilename)
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

    auto getNodeMember = [&](
                             const rapidjson::Value& node,
                             const char* label) -> const rapidjson::Value&
    {
        const rapidjson::Value::ConstMemberIterator iter(node.FindMember(label));
        if (iter == node.MemberEnd()) missingNodeError(thetaFilename, label);
        return iter->value;
    };

    static const char* thetasLabel = "thetas";
    const rapidjson::Value& thetasArray(getNodeMember(document, thetasLabel));
    if (!thetasArray.IsArray()) wrongValueTypeError(thetaFilename, thetasLabel, "array");

    if (thetasArray.Empty())
    {
        std::ostringstream oss;
        oss << "ERROR: No theta values in json indel parameter file '" << thetaFilename << "'\n";
        BOOST_THROW_EXCEPTION(illumina::common::LogicException(oss.str()));
    }

    std::map<unsigned, std::vector<double>> thetasMap;
    for (const auto& thetasByPatternSize : thetasArray.GetArray())
    {
        static const char* repeatPatternSizeLabel = "repeatPatternSize";
        const rapidjson::Value& repeatPatternSizeValue(getNodeMember(thetasByPatternSize, repeatPatternSizeLabel));
        if (!repeatPatternSizeValue.IsUint()) wrongValueTypeError(thetaFilename, repeatPatternSizeLabel, "unsigned");
        const unsigned repeatPatternSize(repeatPatternSizeValue.GetUint());

        static const char* thetaLabel = "theta";
        const rapidjson::Value& thetaArray(getNodeMember(thetasByPatternSize, thetaLabel));
        if (!thetaArray.IsArray()) wrongValueTypeError(thetaFilename, thetaLabel, "array");

        std::vector<double> theta;
        for (const auto& thetaValue : thetaArray.GetArray())
        {
            theta.push_back(thetaValue.GetDouble());
        }
        thetasMap[repeatPatternSize] = theta;
    }
    return thetasMap;
}

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

#pragma once

#include "IndelErrorRateSet.hh"
#include "rapidjson/document.h"
#include "common/RapidJsonHelper.hh"
#include <map>

/// \brief Stores the parameters of the IndelMotifBinomialMixture for a given STR size and repeat length
///
/// Part of the indelErrorModel.json definition
/// Whenever a member variable is added or type is modified, the serialize and deserialize methods must be updated
///
class IndelMotifBinomialMixture
{
public:
    IndelMotifBinomialMixture(const unsigned repeatPatternSize,
                              const unsigned repeatCount,
                              const double indelRate,
                              const double noisyLocusRate):
        _repeatPatternSize(repeatPatternSize)
        , _repeatCount(repeatCount)
        , _indelRate(indelRate)
        , _noisyLocusRate(noisyLocusRate)
    {}
private:
    unsigned _repeatPatternSize = 0;
    unsigned _repeatCount = 0;
    double _indelRate = 0;
    double _noisyLocusRate = 0;
public:
    unsigned getRepeatPatternSize() const
    {
        return _repeatPatternSize;
    }
    unsigned getRepeatCount() const
    {
        return _repeatCount;
    }
    double getIndelRate() const
    {
        return _indelRate;
    }
    double getNoisyLocusRate() const
    {
        return _noisyLocusRate;
    }

    // |brief Serialize object into jsonWriter via rapidJson
    ///
    /// \param[in] writer The writer to serilize into (https://github.com/Tencent/rapidjson/blob/master/example/serialize/serialize.cpp)
    ///
    template <typename Writer>
    void serialize(Writer& writer) const
    {
        writer.String("repeatPatternSize");
        writer.Uint(_repeatPatternSize);
        writer.String("repeatCount");
        writer.Uint(_repeatCount);
        writer.String("indelRate");
        writer.Double(_indelRate);
        writer.String("noisyLocusRate");
        writer.Double(_noisyLocusRate);
    }

    // |brief Deserialize json document to object
    ///
    /// \param[in] root The json document to deserialize
    ///
    static
    IndelMotifBinomialMixture deserialize(const rapidjson::Value& root)
    {
        using namespace illumina::common;

        static const char* indelRateLabel = "indelRate";
        const rapidjson::Value& indelRateValue(RapidJsonHelper::getNodeMember(root, indelRateLabel));
        if (!indelRateValue.IsNumber()) RapidJsonHelper::wrongValueTypeError(indelRateLabel, "number");
        const double indelRate(indelRateValue.GetDouble());

        static const char* noisyLocusRateLabel = "noisyLocusRate";
        const rapidjson::Value& noisyLocusRateValue(RapidJsonHelper::getNodeMember(root, noisyLocusRateLabel));
        if (!noisyLocusRateValue.IsNumber()) RapidJsonHelper::wrongValueTypeError(noisyLocusRateLabel, "number");
        const double noisyLocusRate(noisyLocusRateValue.GetDouble());

        static const char* repeatCountLabel = "repeatCount";
        const rapidjson::Value& repeatCountValue(RapidJsonHelper::getNodeMember(root, repeatCountLabel));
        if (!repeatCountValue.IsUint()) RapidJsonHelper::wrongValueTypeError(repeatCountLabel, "unsigned");
        const unsigned repeatCount(repeatCountValue.GetUint());

        static const char* repeatPatternSizeLabel = "repeatPatternSize";
        const rapidjson::Value& repeatPatternSizeValue(RapidJsonHelper::getNodeMember(root, repeatPatternSizeLabel));
        if (!repeatPatternSizeValue.IsUint()) RapidJsonHelper::wrongValueTypeError(repeatPatternSizeLabel, "unsigned");
        const unsigned repeatPatternSize(repeatPatternSizeValue.GetUint());

        return IndelMotifBinomialMixture(repeatPatternSize, repeatCount, indelRate, noisyLocusRate);
    }
};

/// \brief Stores the list of IndelMotifBinomialMixture objects
///
/// Part of the indelErrorModel.json definition
/// Whenever a member variable is added or type is modified, the serialize and deserialize methods must be updated
///
class IndelErrorModelBinomialMixture
{
public:
    void addMotif(const IndelMotifBinomialMixture& motif)
    {
        _motifs.push_back(motif);
    }
    const std::vector<IndelMotifBinomialMixture>&
    getMotifs() const
    {
        return _motifs;
    }
private:
    std::vector<IndelMotifBinomialMixture> _motifs;
public:
    // |brief Serialize object into jsonWriter via rapidJson
    ///
    /// \param[in] writer The writer to serilize into (https://github.com/Tencent/rapidjson/blob/master/example/serialize/serialize.cpp)
    ///
    template <typename Writer>
    void serialize(Writer& writer) const
    {
        writer.StartArray();

        for (const IndelMotifBinomialMixture& motif:_motifs)
        {
            writer.StartObject();
            motif.serialize(writer);
            writer.EndObject();
        }
        writer.EndArray();
    }

    // |brief Deserialize json document to object
    ///
    /// \param[in] root The json document to deserialize
    ///
    static
    IndelErrorModelBinomialMixture deserialize(const rapidjson::Value& root)
    {
        using namespace illumina::common;

        const rapidjson::Value& motifArray(root);

        IndelErrorModelBinomialMixture indelErrorModelBinomialMixture;
        if (!motifArray.IsArray()) RapidJsonHelper::wrongValueTypeError("motif", "array");
        for (const auto& motifValue : motifArray.GetArray())
        {
            indelErrorModelBinomialMixture.addMotif(IndelMotifBinomialMixture::deserialize(motifValue));
        }

        return indelErrorModelBinomialMixture;
    }

};

/// \brief Stores the contents from a single sample in the indelErrorModel.json in memory
///
/// Part of the indelErrorModel.json definition
/// Whenever a member variable is added or type is modified, the serialize and deserialize methods must be updated
///
class IndelErrorModelJson
{

public:
    /// \brief Initialize json object (used for serialization)
    ///
    /// \param[in] sampleName The sample name (typically the bam file path) used to estimate the model
    ///
    /// \param[in] model The list of motifs estimated from the sample
    ///
    /// \param[in] isStatic Boolean flag indicating whether the model was estimated from the data
    ///
    explicit
    IndelErrorModelJson(const std::string& sampleName, const IndelErrorModelBinomialMixture& model, const bool isStatic);

    // |brief Serialize object into jsonWriter via rapidJson
    ///
    /// \param[in] writer The writer to serilize into (https://github.com/Tencent/rapidjson/blob/master/example/serialize/serialize.cpp)
    ///
    template <typename Writer>
    void serialize(Writer& writer) const
    {
        writer.StartObject();
        writer.String("sampleName");
        writer.String(_sampleName.c_str());
        writer.String("motif");
        _model.serialize(writer);
        writer.String("isStatic");
        writer.Bool(_isStatic);
        writer.EndObject();
    }

    // |brief Deserialize json document to object
    ///
    /// \param[in] root The json document to deserialize
    ///
    static
    IndelErrorModelJson deserialize(const rapidjson::Value& root)
    {
        using namespace illumina::common;

        static const char* sampleNameLabel = "sampleName";
        const rapidjson::Value& sampleNameValue(RapidJsonHelper::getNodeMember(root, sampleNameLabel));
        if (!sampleNameValue.IsString()) RapidJsonHelper::wrongValueTypeError(sampleNameLabel, "string");
        const std::string sampleName(sampleNameValue.GetString());

        static const char* motifLabel = "motif";
        const rapidjson::Value& motifArray(RapidJsonHelper::getNodeMember(root, motifLabel));
        IndelErrorModelBinomialMixture model(IndelErrorModelBinomialMixture::deserialize(motifArray));

        static const char* isStaticLabel = "isStatic";
        const rapidjson::Value& isStaticValue(RapidJsonHelper::getNodeMember(root, isStaticLabel));
        if (!isStaticValue.IsBool()) RapidJsonHelper::wrongValueTypeError(isStaticLabel, "bool");
        const bool isStatic(isStaticValue.GetBool());

        return IndelErrorModelJson(sampleName, model, isStatic);
    }

    const std::string&
    getSampleName() const
    {
        return _sampleName;
    }

    void
    setSampleName(const std::string& sampleName)
    {
        _sampleName = sampleName;
    }

    const IndelErrorModelBinomialMixture&
    getBinomialMixtureModel() const
    {
        return _model;
    }

private:
    std::string _sampleName;
    IndelErrorModelBinomialMixture _model;
    bool _isStatic;
};


/// \brief Stores the contents from the indelErrorModel.json in memory
///
/// Part of the indelErrorModel.json definition
/// Whenever a member variable is added or type is modified, the serialize and deserialize methods must be updated
///
class IndelErrorModelsJson
{
public:
    IndelErrorModelsJson() {}

    // |brief Serialize object into jsonWriter via rapidJson
    ///
    /// \param[in] writer The writer to serilize into (https://github.com/Tencent/rapidjson/blob/master/example/serialize/serialize.cpp)
    ///
    template <typename Writer>
    void serialize(Writer& writer) const
    {
        writer.StartObject();
        writer.String("sample");
        writer.StartArray();
        for (const auto& model:_models)
        {
            model.serialize(writer);
        }
        writer.EndArray();
        writer.EndObject();
    }

    // |brief Deserialize json document to object
    ///
    /// \param[in] root The json document to deserialize
    ///
    static
    IndelErrorModelsJson deserialize(const rapidjson::Value& root)
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

    void addModel(const IndelErrorModelJson& model)
    {
        _models.push_back(model);
    }

    const std::vector<IndelErrorModelJson>& getIndelErrorModels() const
    {
        return _models;
    }

    IndelErrorModelJson& getIndelErrorModel(const unsigned modelIndex)
    {
        assert(modelIndex < _models.size());
        return _models[modelIndex];
    }

private:
    std::vector<IndelErrorModelJson> _models;
};

/// \brief Handles all serialization/deserialization to/from json files for the IndelErrorModel
///
class IndelErrorModelParser
{
public:
    /// \brief Deserializes the indel error rate values for each sample
    ///
    /// \param[in] modelFilename The json filename to deserialize
    /// \param[in] modelMap The map to import into
    ///
    static void
    importIndelErrorModelJsonFile(
        const std::string& modelFilename,
        std::map<std::string, IndelErrorRateSet>& modelMap);
    /// \brief Deserializes the indel error rate values for each sample
    ///
    /// \param[in] modelFilename The json filename to deserialize
    /// \param[in] indelErrorModelsJson The object to import into
    ///
    static void
    importIndelErrorModelJsonFile(
        const std::string& modelFilename,
        IndelErrorModelsJson& indelErrorModelsJson);
    /// \brief Deserializes the theta values for each repeat pattern size
    ///
    /// \param[in] thetaFilename The json filename to deserialize
    /// \param[in] thetasMap The map to import into
    ///
    static void
    importThetaJsonFile(
        const std::string& thetaFilename,
        std::map<unsigned, std::vector<double>>& thetasMap);

    /// \brief Deserializes multiple json files and populates the IndelErrorRateSet object for each sample
    ///
    /// \param[in] modelFilenames The json filenames to deserialize
    ///
    static std::map<std::string, IndelErrorRateSet>
    generateIndelErrorRateSetMap(
        const std::vector<std::string>& modelFilenames);
};

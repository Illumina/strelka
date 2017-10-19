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

#pragma once

#include "IndelErrorRateSet.hh"
#include "rapidjson/document.h"
#include "common/RapidJsonHelper.hh"
#include <map>


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
    static
    IndelMotifBinomialMixture deserialize(const rapidjson::Value& root)
    {
        using namespace illumina::common;

        static const char* indelRateLabel = "indelRate";
        const rapidjson::Value& indelRateValue(RapidJsonHelper::getNodeMember(root, indelRateLabel));
        const double indelRate(indelRateValue.GetDouble());

        static const char* noisyLocusRateLabel = "noisyLocusRate";
        const rapidjson::Value& noisyLocusRateValue(RapidJsonHelper::getNodeMember(root, noisyLocusRateLabel));
        const double noisyLocusRate(noisyLocusRateValue.GetDouble());

        static const char* repeatCountLabel = "repeatCount";
        const rapidjson::Value& repeatCountValue(RapidJsonHelper::getNodeMember(root, repeatCountLabel));
        const unsigned repeatCount(repeatCountValue.GetUint());

        static const char* repeatPatternSizeLabel = "repeatPatternSize";
        const rapidjson::Value& repeatPatternSizeValue(RapidJsonHelper::getNodeMember(root, repeatPatternSizeLabel));
        const unsigned repeatPatternSize(repeatPatternSizeValue.GetUint());

        return IndelMotifBinomialMixture(repeatPatternSize, repeatCount, indelRate, noisyLocusRate);
    }
};

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
    static
    IndelErrorModelBinomialMixture deserialize(const rapidjson::Value& root)
    {
        using namespace illumina::common;

        //static const char* motifLabel = "motif";
        const rapidjson::Value& motifArray(root);

        IndelErrorModelBinomialMixture indelErrorModelBinomialMixture;
        for (const auto& motifValue : motifArray.GetArray())
        {
            indelErrorModelBinomialMixture.addMotif(IndelMotifBinomialMixture::deserialize(motifValue));
        }

        return indelErrorModelBinomialMixture;
    }

};


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

    static IndelErrorModelJson deserialize(const rapidjson::Value& root);

    const std::string&
    getSampleName() const
    {
        return _sampleName;
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

class IndelErrorModelsJson
{
public:
    IndelErrorModelsJson() {}
    static IndelErrorModelsJson deserialize(const rapidjson::Value& root);
    void addModel(const IndelErrorModelJson& model)
    {
        _models.push_back(model);
    }
    const std::vector<IndelErrorModelJson>& getIndelErrorModels() const
    {
        return _models;
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

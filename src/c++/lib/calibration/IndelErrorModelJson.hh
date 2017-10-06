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
#include <map>
#include "rapidjson/document.h"
#include "IndelErrorRateSet.hh"

class IndelMotifBinomialMixture
{
public:
    IndelMotifBinomialMixture(){}
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
};

class IndelErrorModelBinomialMixture
{
public:
    IndelErrorModelBinomialMixture(){};
    void addMotif(const IndelMotifBinomialMixture &motif)
    {
        _motifs.push_back(motif);
    }
private:
    std::vector<IndelMotifBinomialMixture> _motifs;
public:
    template <typename Writer>
    void serialize(Writer& writer) const
    {
        writer.StartArray();

        for(const IndelMotifBinomialMixture &motif:_motifs)
        {
            writer.StartObject();
            motif.serialize(writer);
            writer.EndObject();
        }
        writer.EndArray();
    }
};


/// \brief Handles all serialization/deserialization to/from json files for the IndelErrorModel
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


    /// \brief Deserializes multiple json files and populates the IndelErrorRateSet object for each sample
    ///
    /// \param[in] modelFilenames The json filenames to deserialize
    ///
    static std::map<std::string, IndelErrorRateSet>
    deserializeIndelErrorModels(
        const std::vector<std::string>& modelFilenames);

    /// \brief Deserializes the theta values for each repeat pattern size
    ///
    /// \param[in] filename The json filename to deserialize
    ///
    static std::map<unsigned, std::vector<double> >
    deserializeTheta(
        const std::string& filename);


    std::string
    getSampleName() const
    {
        return _sampleName;
    }

private:
    std::string _sampleName;
    IndelErrorModelBinomialMixture _model;
    bool _isStatic;

};





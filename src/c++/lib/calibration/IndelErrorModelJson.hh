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

#include "json/json.h"
#include "IndelErrorRateSet.hh"

class IndelMotifBinomialMixture
{
public:
    unsigned repeatPatternSize = 0;
    unsigned repeatCount = 0;
    double indelRate = 0;
    double noisyLocusRate = 0;
};

class IndelModelBinomialMixture
{
public:
    std::vector<IndelMotifBinomialMixture> motifs;
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
    explicit
    IndelErrorModelJson(const std::string& sampleName);

    /// \brief Adds the motif with the given parameters to the json object
    ///
    /// \param[in] repeatPatternSize The length of the repeat pattern
    ///
    /// \param[in] repeatCount Number of repetitions for the given repeatPatternSize
    ///
    /// \param[in] indelRate The estimated indel error rate
    ///
    /// \param[in] noisyLocusRate The probability that a locus is in a noisy state
    ///
    void addMotif(
        unsigned repeatPatternSize,
        unsigned repeatCount,
        double indelRate,
        double noisyLocusRate);

    // TODO: This isn't really serialization. Serialize the IndelErrorRateSet instead
    /// \brief Serializes the model and writes it out to a json file
    ///
    /// \param[in] sampleName The sample name (typically the bam file path) used to estimate the model
    ///
    /// \param[in] motifsNode The json value to write out
    ///
    /// \param[in] isStatic Flag describing whether the model params were estimated from the specific sample or taken from the static model
    ///
    /// \param[in] filename The name of the json file to write to
    ///
    static void
    serializeIndelErrorModel(
        const std::string& sampleName,
        const Json::Value& motifsNode,
        const bool isStatic,
        const std::string& filename);

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

    Json::Value
    generateMotifsNode() const;

    std::string
    getSampleName() const
    {
        return _sampleName;
    }

public:
    IndelModelBinomialMixture model;

private:
    std::string _sampleName;


};





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


class IndelErrorModelJson
{

public:
    explicit
    IndelErrorModelJson(const std::string& sampleName);

    void addMotif(
            unsigned repeatPatternSize,
            unsigned repeatCount,
            double indelRate,
            double noisyLocusRate);

    void exportIndelErrorModelToJsonFile(
            const std::string& filename) const;

    static void
    writeIndelErrorModelJsonFile(
            const std::string& sampleName,
            const Json::Value& motifsNode,
            const std::string& filename);

    static std::map<std::string, IndelErrorRateSet>
    deserializeIndelModels(
            const std::vector<std::string>& modelFilenames);

    static std::map<unsigned, std::vector<double> >
    importTheta(
            const std::string& filename);

public:
    IndelModelBinomialMixture model;

private:
    std::string _sampleName;

    Json::Value
    generateMotifsNode() const;
};





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

#include "calibration/IndelErrorModel.hh"
#include "json/json.h"
#include "errorAnalysis/SequenceErrorCounts.hh"
class IndelModelProduction
{
public:
    IndelModelProduction(
        const SequenceErrorCounts& counts,
        const std::string& thetaFilename,
        const std::string& outputFilename);
    void estimateIndelErrorRates();

    void exportModel() const;

    void exportModelUsingInputJson(
        const std::string& jsonFilename) const;

    bool checkEstimatedModel() const;

private:
    std::map<unsigned, std::vector<double> >
    importTheta(
        const std::string& filename);
    bool isValidErrorRate(
        const double indelErrorRate) const;

private:
    bool _isEstimated = false;
    SequenceErrorCounts _counts;
    std::string _outputFilename;
    std::map<unsigned, std::vector<double> > _thetas;
    std::vector<AdaptiveIndelErrorModel> _adaptiveIndelErrorModels;
    std::vector<unsigned> _repeatPatterns = {1, 2};
    std::vector<unsigned> _maxRepeatCounts = {16, 9};
    AdaptiveIndelErrorModelLogParams _nonSTRModelParams;
    const double _maxErrorRate = 1.0;
};

// move these to a more appropriate place later
// TODO: these classes can be automatically serialized with cereal
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

class IndelModelJson
{
public:
    IndelModelBinomialMixture model;

    explicit
    IndelModelJson(const std::string& sampleName);

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

private:
    std::string _sampleName;

    Json::Value
    generateMotifsNode() const;


};

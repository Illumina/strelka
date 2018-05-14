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

#include "calibration/IndelErrorModelJson.hh"
#include "calibration/IndelErrorModel.hh"
#include "errorAnalysis/SequenceAlleleCounts.hh"


class IndelModelProduction
{
public:
    IndelModelProduction(
        const SequenceAlleleCounts& counts,
        const std::string& thetaFilename,
        const std::string& outputFilename);

    void estimateIndelErrorRates();

    IndelErrorModelJson generateIndelErrorModelJson() const;

    void exportIndelErrorModelJson() const;

    void exportModelUsingInputJson(
        const std::string& modelFilename) const;

    static void
    writeIndelErrorModelJson(
        const IndelErrorModelsJson& indelErrorModelsJson,
        const std::string& outputFileName);

    bool checkEstimatedModel() const;

private:
    bool isValidErrorRate(
        const double indelErrorRate) const;

private:
    bool _isEstimated = false;
    bool _isEstimationAcceptable = true;
    SequenceAlleleCounts _counts;
    std::string _outputFilename;
    std::map<unsigned, std::vector<double> > _thetas;
    std::vector<AdaptiveIndelErrorModel> _adaptiveIndelErrorModels;
    std::vector<unsigned> _repeatPatterns = {1, 2};
    std::vector<unsigned> _maxRepeatCounts = {16, 9};
    AdaptiveIndelErrorModelLogParams _nonSTRModelParams;
    const double _maxErrorRate = .3;
    const double _minErrorRate = .0;
};

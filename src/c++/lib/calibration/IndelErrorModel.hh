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
#include "IndelErrorModelMetadata.hh"
#include "IndelErrorRateSet.hh"

#include "starling_common/AlleleReportInfo.hh"


/// \brief Organizes indel error rate information.
///
struct IndelErrorModel
{
    /// \brief Initialize indel error model to either a precomputed static model (if \p modelFilename is empty),
    /// or from a json parameter file otherwise.
    ///
    IndelErrorModel(
        const std::string& modelName,
        const std::string& modelFilename);

    /// \brief Retrieve indel error rates for a specific indel type.
    ///
    /// \param[in] isCandidateRates If true, retrieve rates to be used for indel candidate testing.
    void
    getIndelErrorRate(
        const IndelKey& indelKey,
        const AlleleReportInfo& indelReportInfo,
        double& refToIndelErrorProb,
        double& indelToRefErrorProb,
        const bool isCandidateRates = false) const;

    void
    deserializeLegacyIndelModels(
        const std::string& modelName,
        const std::string& modelFilename,
        const Json::Value& root);

    void
    deserializeIndelModels
    (const std::string& modelFilename,
     const Json::Value& root);

private:
    IndelErrorModelMetadata _meta;
    IndelErrorRateSet _errorRates;

    /// error rates used for candidate indel selection only
    IndelErrorRateSet _candidateErrorRates;
};

class AdaptiveIndelErrorModelLogParams
{
public:
    double logErrorRate = -std::numeric_limits<double>::infinity();
    double logNoisyLocusRate = -std::numeric_limits<double>::infinity();
};

// TODO: This class will be useful when we put in the production estimator
class AdaptiveIndelErrorModel
{
public:
    AdaptiveIndelErrorModel(
        unsigned repeatPatternSize,
        unsigned highRepeatCount,
        const AdaptiveIndelErrorModelLogParams& lowLogParams,
        const AdaptiveIndelErrorModelLogParams& highLogParams);
private:
    unsigned _repeatPatternSize = 0;
    unsigned _lowRepeatCount = lowRepeatCount;
    unsigned _highRepeatCount = 0;

    AdaptiveIndelErrorModelLogParams _lowLogParams;
    AdaptiveIndelErrorModelLogParams _highLogParams;

public:
    unsigned
    repeatPatternSize() const
    {
        return _repeatPatternSize;
    }
    unsigned
    highRepeatCount() const
    {
        return _highRepeatCount;
    }

    double
    errorRate(
        const unsigned repeatCount) const;
    double
    noisyLocusRate(
        const unsigned repeatCount) const;


    /// Perform a linear fit from 2 known points and return y corresponding to x
    ///
    /// \param x the point on the linear curve of interest
    /// \param x1 the first known x position
    /// \param y1 the y position corresponding to x1
    /// \param x2 the 2nd known x position
    /// \param y2 the y position corresponding to x2
    static double
    linearFit(
        const double x,
        const double x1,
        const double y1,
        const double x2,
        const double y2);

    static unsigned lowRepeatCount; // it should be safe to fix this to 2
};

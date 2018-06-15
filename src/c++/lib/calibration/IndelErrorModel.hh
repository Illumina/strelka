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

#include "starling_common/AlleleReportInfo.hh"


/// \brief Organizes indel error rate information.
///
struct IndelErrorModel
{
    /// \brief Initialize indel error model to either a precomputed static model (if \p modelFilenames is empty),
    /// or from a json parameter file otherwise.
    ///
    /// \param[in] alignmentFilenames Name and indexed order of alignment files, which can be used to sync sample
    ///                           index values in the indel error model file
    /// \param[in] modelName Name of selected static indel error model to use, ignored if \p modelFilenames is non-empty
    /// \param[in] modelFilenames Indel error model params in json format, one json file per sample.
    IndelErrorModel(
        const std::vector<std::string>& alignmentFilenames,
        const std::string& modelName,
        const std::vector<std::string>& modelFilenames);

    /// \brief Retrieve indel error rates for a specific indel type.
    ///
    /// The function currently uses the sample index as a proxy for a set of read groups with a shared sample-prep/
    /// sequencing error noise profile. Over time sampleIndex should be replaced with a more correct type of
    /// 'sequencingAssayGroupIndex'
    ///
    /// \param[in] sampleIndex Index of the sample for which indel error rates should be returned
    /// \param[in] isCandidateRates If true, retrieve rates to be used for indel candidate testing
    void
    getIndelErrorRate(
        const unsigned sampleIndex,
        const IndelKey& indelKey,
        const AlleleReportInfo& indelReportInfo,
        double& refToIndelErrorProb,
        double& indelToRefErrorProb,
        const bool isCandidateRates = false) const;

private:
    void
    checkSampleIndex(
        const unsigned sampleIndex) const;

    IndelErrorRateSet&
    getSampleSpecificIndelErrorRates(
        const unsigned sampleIndex = 0)
    {
        checkSampleIndex(sampleIndex);
        const unsigned sampleIndexUsed(_isUseSampleSpecificErrorRates ? sampleIndex : 0);
        return _sampleErrorRates[sampleIndexUsed];
    }

    const IndelErrorRateSet&
    getSampleSpecificIndelErrorRates(
        const unsigned sampleIndex) const
    {
        checkSampleIndex(sampleIndex);
        const unsigned sampleIndexUsed(_isUseSampleSpecificErrorRates ? sampleIndex : 0);
        return _sampleErrorRates[sampleIndexUsed];
    }

    bool _isUseSampleSpecificErrorRates = false;

    /// \brief Standard indel error rates for each sample
    std::vector<IndelErrorRateSet> _sampleErrorRates;

    /// \brief Error rates used for candidate indel selection only
    IndelErrorRateSet _candidateErrorRates;

    /// \brief Track total number of expected samples in support of
    /// sample-specific error rates
    unsigned _sampleCount;
};


class AdaptiveIndelErrorModelLogParams
{
public:
    double logErrorRate = -std::numeric_limits<double>::infinity();
    double logNoisyLocusRate = -std::numeric_limits<double>::infinity();
    bool paramsAcceptable = true;
};


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

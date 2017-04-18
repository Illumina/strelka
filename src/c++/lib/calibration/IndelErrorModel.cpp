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
/*
 *      Author: Morten Kallberg
 */

#include "IndelErrorModel.hh"
#include "blt_util/log.hh"
#include "common/Exceptions.hh"

#include "json/json.h"

#include <cmath>
#include <cassert>

#include <fstream>

/// simple log-linear error ramp as a function of hpol length - default error model used in NS5/v2.7.x release series
static
IndelErrorRateSet
getLogLinearIndelErrorModel()
{
    static const double logLowErrorRate(std::log(5e-5));
    static const double logHighErrorRate(std::log(3e-4));

    // this is the zero-indexed endpoint of the ramp, so we hit the
    // constant high error rate at an hpol length of repeatCountSwitchPoint+1
    static const unsigned repeatCountSwitchPoint(15);

    IndelErrorRateSet rates;

    // model covers homopolymers only:
    static const unsigned repeatingPatternSize(1);

    for (unsigned patternRepeatCount=1; patternRepeatCount <= (repeatCountSwitchPoint+1); ++patternRepeatCount)
    {
        const double highErrorFrac(std::min((patternRepeatCount-1),repeatCountSwitchPoint)/static_cast<double>(repeatCountSwitchPoint));
        const double logErrorRate((1.-highErrorFrac)*logLowErrorRate + highErrorFrac*logHighErrorRate);
        const double errorRate(std::exp(repeatingPatternSize==1 ? logErrorRate : logLowErrorRate));

        rates.addRate(repeatingPatternSize, patternRepeatCount, errorRate, errorRate);
    }
    return rates;
}



/// Reads in model parameter matrix with entries as error-pair [ins_error,del_error]
/// in the following format:
/// unit length 1: [[[del_hpol1,ins_hpol1],[del_hpol2,ins_hpol2],...,[del_hpol_m,ins_hpol_m]],
/// unit length 2:  [[del_dinuc1,ins_dinuc1],[del_dinuc2,ins_dinuc2],...,[del_dinuc_m,ins_dinuc_m]],
///  ....
/// unit length N:  [[del_repeatN,ins_repeatN],[del_repeatN2,ins_repeatN2],...,]]]
static
IndelErrorRateSet
deserializeRateSet(
    const Json::Value& root)
{
    const unsigned maxRepeatingPatternSize = root["MaxMotifLength"].asInt();
    const unsigned maxTractLength = root["MaxTractLength"].asInt();
    const Json::Value& models = root["Model"];

    assert((models.size()==maxRepeatingPatternSize) && "Unexpected repeating pattern size in indel model");

    IndelErrorRateSet rates;

    for (unsigned repeatingPatternSize(1); repeatingPatternSize<=maxRepeatingPatternSize; ++repeatingPatternSize)
    {
        const auto& pattern(models[repeatingPatternSize-1]);

        assert((pattern.size()<=maxTractLength) && "Unexpected tract length in indel model");

        for (unsigned tractLength=1; tractLength <= pattern.size(); ++tractLength)
        {
            if ((tractLength % repeatingPatternSize)!=0) continue;
            const unsigned patternRepeatCount(tractLength/repeatingPatternSize);
            const auto& cell(pattern[tractLength-1]);

            const double delete_error_prob(cell[0].asDouble());
            const double insert_error_prob(cell[1].asDouble());
            rates.addRate(repeatingPatternSize, patternRepeatCount, insert_error_prob, delete_error_prob);
        }
    }

    return rates;
}



IndelErrorModel::
IndelErrorModel(
    const std::string& modelName,
    const std::string& modelFilename)
{
    if (modelFilename.empty())
    {
        if (modelName == "logLinear")
        {
            _errorRates = getLogLinearIndelErrorModel();
        }
        else
        {
            using namespace illumina::common;

            std::ostringstream oss;
            oss << "ERROR: unrecognized indel error model name: '" << modelName << "'\n";
            BOOST_THROW_EXCEPTION(LogicException(oss.str()));
        }
    }
    else
    {
        Json::Value root;
        std::ifstream file(modelFilename , std::ifstream::binary);
        file >> root;

        Json::Value models = root["IndelModels"];
        if (! models.isNull())
        {
            for (const auto& modelValue : models)
            {
                _meta.deserialize(modelValue);
                if (_meta.name != modelName) continue;
                _errorRates = deserializeRateSet(modelValue);
                break;
            }
        }

        if (_meta.name != modelName)
        {
            using namespace illumina::common;

            std::ostringstream oss;
            oss << "ERROR: unrecognized indel error model name: '" << modelName << "' in model file '" << modelFilename << "'\n";
            BOOST_THROW_EXCEPTION(LogicException(oss.str()));
        }
    }

    _errorRates.finalizeRates();
}




void
IndelErrorModel::
getIndelErrorRate(
    const IndelKey& indelKey,
    const AlleleReportInfo& indelReportInfo,
    double& refToIndelErrorProb,
    double& indelToRefErrorProb) const
{
    using namespace IndelErrorRateType;

    const index_t indelType(getRateType(indelKey));
    // determine simple case
    const bool isSimpleIndel(indelType==INSERT || indelType==DELETE);

    if (! isSimpleIndel)
    {
        // complex indels use baseline indel error rates
        /// TODO - provide estimates for complex indels
        const double baselineInsertionErrorRate(_errorRates.getRate(1,1,INSERT));
        const double baselineDeletionErrorRate(_errorRates.getRate(1,1,DELETE));

        refToIndelErrorProb=std::max(baselineInsertionErrorRate,baselineDeletionErrorRate);
        indelToRefErrorProb=refToIndelErrorProb;
        return;
    }
    else
    {
        // determine the repeat pattern size and count:
        static const unsigned one(1);
        const unsigned repeatingPatternSize = std::max(indelReportInfo.repeat_unit_length, one);
        const unsigned refPatternRepeatCount = std::max(indelReportInfo.ref_repeat_count, one);
        const unsigned indelPatternRepeatCount = std::max(indelReportInfo.indel_repeat_count, one);

        //const int indelPatternRepeatSize(std::abs(static_cast<int>(iri.ref_repeat_count)-static_cast<int>(iri.indel_repeat_count));

        const index_t reverseIndelType((indelType == DELETE) ? INSERT : DELETE);

        refToIndelErrorProb = _errorRates.getRate(repeatingPatternSize, refPatternRepeatCount, indelType);
        indelToRefErrorProb = _errorRates.getRate(repeatingPatternSize, indelPatternRepeatCount, reverseIndelType);
    }
}

// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2016 Illumina, Inc.
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


// Calculate p(error) of
//    CASE: del
//    FIT pars: [  1.49133831e-03   1.03348683e+01   1.13646811e+00   1.18488282e-05]
//    Function prob(error)=0.00149133830825/ (1 + exp((10.3348683003-x)/1.13646810558))+1.18488281756e-05
//    --------------------
//    CASE: ins
//    FIT pars: [  1.09573511e-03   9.82226042e+00   1.03579658e+00   8.31843836e-06]
//    Function prob(error)=0.00109573511176/ (1 + exp((9.82226041538-x)/1.03579658224))+8.31843836296e-06
//    --------------------
static
IndelErrorRateSet
getNewIndelErrorModel()
{
    static const unsigned maxPatternRepeatCount = 40;

    static const double insert_A(1.49133831e-03);
    static const double insert_B(1.03348683e+01);
    static const double insert_C(1.13646811e+00);
    static const double insert_D(1.18488282e-05);

    static const double delete_A(1.09573511e-03);
    static const double delete_B(9.82226042e+00);
    static const double delete_C(1.03579658e+00);
    static const double delete_D(8.31843836e-06);

    IndelErrorRateSet rates;

    // model covers homopolymers only:
    static const unsigned repeatingPatternSize(1);

    for (unsigned patternRepeatCount=1; patternRepeatCount <= maxPatternRepeatCount; ++patternRepeatCount)
    {
        const double insert_g(insert_A/ (1 + std::exp((insert_B-patternRepeatCount)/insert_C))+insert_D);
        const double insert_error_prob(1.-std::exp(-insert_g/patternRepeatCount));
        const double delete_g(delete_A/ (1 + std::exp((delete_B-patternRepeatCount)/delete_C))+delete_D);
        const double delete_error_prob(1.-std::exp(-delete_g/patternRepeatCount));

        rates.addRate(repeatingPatternSize, patternRepeatCount, insert_error_prob, delete_error_prob);
    }
    return rates;
}



static
IndelErrorRateSet
getOldIndelErrorModel()
{
    static const unsigned maxPatternRepeatCount = 40;

    static const double insert_A(5.03824e-7);
    static const double insert_B(3.30572e-10);
    static const double insert_C(6.99777);

    static const double delete_hpol1_err(5.00057e-5);
    static const double delete_A(1.09814e-5);
    static const double delete_B(5.19742e-10);
    static const double delete_C(6.99256);

    IndelErrorRateSet rates;

    // model covers homopolymers only:
    static const unsigned repeatingPatternSize(1);

    for (unsigned patternRepeatCount=1; patternRepeatCount <= maxPatternRepeatCount; ++patternRepeatCount)
    {
        const double insert_g(insert_A*patternRepeatCount+insert_B*std::pow(patternRepeatCount,insert_C));
        const double insert_error_prob(1.-std::exp(-insert_g));

        double delete_g(delete_hpol1_err);
        if (patternRepeatCount>1)
        {
            delete_g = delete_A*patternRepeatCount+delete_B*std::pow(patternRepeatCount,delete_C);
        }
        const double delete_error_prob(1.-std::exp(-delete_g));

        rates.addRate(repeatingPatternSize, patternRepeatCount, insert_error_prob, delete_error_prob);
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
        if     (modelName == "old")
        {
            _errorRates = getOldIndelErrorModel();
        }
        else if(modelName == "new")
        {
            _errorRates = getNewIndelErrorModel();
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

        if(_meta.name != modelName)
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
    const starling_indel_report_info& iri,
    double& refToIndelErrorProb,
    double& indelToRefErrorProb) const
{
    // determine simple case
    const bool isSimpleIndel(iri.it==INDEL::INSERT || iri.it==INDEL::DELETE);

    if (! isSimpleIndel)
    {
        // complex indels use baseline indel error rates
        /// TODO - provide estimates for complex indels
        const double baselineInsertionErrorRate(_errorRates.getRate(1,1,INDEL::INSERT));
        const double baselineDeletionErrorRate(_errorRates.getRate(1,1,INDEL::DELETE));

        refToIndelErrorProb=std::max(baselineInsertionErrorRate,baselineDeletionErrorRate);
        indelToRefErrorProb=refToIndelErrorProb;
        return;
    }

    assert(iri.it == INDEL::INSERT || iri.it == INDEL::DELETE);

    // determine the repeat pattern size and count:
    static const unsigned one(1);
    const unsigned repeatingPatternSize    = std::max(iri.repeat_unit_length,one);
    const unsigned refPatternRepeatCount   = std::max(iri.ref_repeat_count,one);
    const unsigned indelPatternRepeatCount = std::max(iri.indel_repeat_count,one);

    //const int indelPatternRepeatSize(std::abs(static_cast<int>(iri.ref_repeat_count)-static_cast<int>(iri.indel_repeat_count));

    const INDEL::index_t reverse_it(iri.it==INDEL::DELETE ? INDEL::INSERT : INDEL::DELETE);

    refToIndelErrorProb=_errorRates.getRate(repeatingPatternSize, refPatternRepeatCount, iri.it);
    indelToRefErrorProb=_errorRates.getRate(repeatingPatternSize, indelPatternRepeatCount, reverse_it);
}

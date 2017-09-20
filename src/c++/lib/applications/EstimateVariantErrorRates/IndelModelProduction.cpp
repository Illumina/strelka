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

#include <iomanip>
#include <iostream>
#include <fstream>

#include "common/Exceptions.hh"
#include "blt_util/log.hh"
#include "blt_util/math_util.hh"
#include "blt_util/prob_util.hh"

#include "IndelModelProduction.hh"
#include "calibration/IndelErrorModelJson.hh"

//#define CODEMIN_DEBUG
#define CODEMIN_USE_BOOST
#include "minimize_conj_direction.h"

namespace MIN_PARAMS3
{
enum index_t
{
    LN_INSERT_ERROR_RATE,
    LN_DELETE_ERROR_RATE,
    LN_NOISY_LOCUS_RATE,
    LN_THETA,
    SIZE
};
}

static
double
getObsLogLhood(
    const double logHomPrior,
    const double logHetPrior,
    const double logAltHetPrior,
    const double logNoIndelPrior,
    const double logInsertErrorRate,
    const double logDeleteErrorRate,
    const double logNoIndelRefRate,
    const ExportedIndelObservations& obs)
{
    static const double log0(-std::numeric_limits<double>::infinity());

    static const double homAltRate(0.99);
    static const double hetAltRate(0.5);

    static const double logHomAltRate(std::log(homAltRate));
    static const double logHomRefRate(std::log(1.-homAltRate));
    static const double logHetRate(std::log(hetAltRate));

    // get lhood of homref GT:
    double noindel(log0);
    {
        unsigned totalInsertObservations(0);
        for (unsigned altIndex(INDEL_SIGNAL_TYPE::INSERT_1); altIndex<INDEL_SIGNAL_TYPE::DELETE_1; ++altIndex)
        {
            totalInsertObservations += obs.altObservations[altIndex];
        }

        unsigned totalDeleteObservations(0);
        for (unsigned altIndex(INDEL_SIGNAL_TYPE::DELETE_1); altIndex<INDEL_SIGNAL_TYPE::SIZE; ++altIndex)
        {
            totalDeleteObservations += obs.altObservations[altIndex];
        }

        noindel = (
                      logInsertErrorRate*totalInsertObservations +
                      logDeleteErrorRate*totalDeleteObservations +
                      logNoIndelRefRate*obs.refObservations);
    }

    unsigned maxIndex(0);
    for (unsigned altIndex(1); altIndex<INDEL_SIGNAL_TYPE::SIZE; ++altIndex)
    {
        if (obs.altObservations[altIndex] > obs.altObservations[maxIndex]) maxIndex = altIndex;
    }

    // get lhood of het and hom GT:
    double het(log0);
    double hom(log0);
    {
        // approximate that the most frequent observations is the only potential variant allele:

        unsigned remainingInsertObservations(0);
        for (unsigned altIndex(INDEL_SIGNAL_TYPE::INSERT_1); altIndex<INDEL_SIGNAL_TYPE::DELETE_1; ++altIndex)
        {
            if (altIndex==maxIndex) continue;
            remainingInsertObservations += obs.altObservations[altIndex];
        }

        unsigned remainingDeleteObservations(0);
        for (unsigned altIndex(INDEL_SIGNAL_TYPE::DELETE_1); altIndex<INDEL_SIGNAL_TYPE::SIZE; ++altIndex)
        {
            if (altIndex==maxIndex) continue;
            remainingDeleteObservations += obs.altObservations[altIndex];
        }

        // compute lhood of het/hom states given that maxIndex is the variant allele:
        het =(logHetRate*(obs.refObservations+obs.altObservations[maxIndex]) +
              logInsertErrorRate*remainingInsertObservations +
              logDeleteErrorRate*remainingDeleteObservations);

        hom = (logHomAltRate*obs.altObservations[maxIndex] +
               logHomRefRate*obs.refObservations +
               logInsertErrorRate*remainingInsertObservations +
               logDeleteErrorRate*remainingDeleteObservations);
    }

    // get lhood of althet GT:
    double althet(log0);
    {
        // approximate that the two most frequent observations are the only potential variant alleles:
        assert(INDEL_SIGNAL_TYPE::SIZE>1);
        unsigned maxIndex2(maxIndex==0 ? 1 : 0);
        for (unsigned altIndex(maxIndex2+1); altIndex<INDEL_SIGNAL_TYPE::SIZE; ++altIndex)
        {
            if (altIndex==maxIndex) continue;
            if (obs.altObservations[altIndex] > obs.altObservations[maxIndex2]) maxIndex2 = altIndex;
        }

        unsigned remainingInsertObservations(0);
        for (unsigned altIndex(INDEL_SIGNAL_TYPE::INSERT_1); altIndex<INDEL_SIGNAL_TYPE::DELETE_1; ++altIndex)
        {
            if (altIndex==maxIndex) continue;
            if (altIndex==maxIndex2) continue;
            remainingInsertObservations += obs.altObservations[altIndex];
        }

        unsigned remainingDeleteObservations(0);
        for (unsigned altIndex(INDEL_SIGNAL_TYPE::DELETE_1); altIndex<INDEL_SIGNAL_TYPE::SIZE; ++altIndex)
        {
            if (altIndex==maxIndex) continue;
            if (altIndex==maxIndex2) continue;
            remainingDeleteObservations += obs.altObservations[altIndex];
        }

        // compute lhood of het/hom states given that maxIndex is the variant allele:
        althet =(logHetRate*(obs.altObservations[maxIndex]+obs.altObservations[maxIndex2]) +
                 logHomRefRate*obs.refObservations +
                 logInsertErrorRate*remainingInsertObservations +
                 logDeleteErrorRate*remainingDeleteObservations);
    }

    return log_sum( log_sum(logHomPrior+hom,logHetPrior+het), log_sum(logNoIndelPrior+noindel,logAltHetPrior+althet) );
}

static
double
contextLogLhood(
    const std::vector<ExportedIndelObservations>& observations,
    const double logInsertErrorRate,
    const double logDeleteErrorRate,
    const double logNoisyLocusRate,
    const double logTheta)
{
#ifdef DEBUG_MODEL3
    log_os << "MODEL3: loghood input:"
           << " insert: " << std::exp(logInsertErrorRate)
           << " delete: " << std::exp(logDeleteErrorRate)
           << " noise: " << std::exp(logNoisyLocusRate)
           << " theta: " << std::exp(logTheta)
           << "\n";
#endif

    static const double log2(std::log(2));
    const double logHomPrior(logTheta-log2);
    const double logHetPrior(logTheta);
    const double logAltHetPrior(logTheta*2);
    const double theta(std::exp(logTheta));
    const double logNoIndelPrior(std::log(1-(theta*3./2.+(theta*theta))));

    const double logNoIndelRefRate(std::log(1-std::exp(logInsertErrorRate)-std::exp(logDeleteErrorRate)));

    const double logCleanLocusRate(std::log(1-std::exp(logNoisyLocusRate)));

    static const double cleanLocusIndelRate(1e-8);
    static const double logCleanLocusIndelRate(std::log(cleanLocusIndelRate));
    static const double logCleanLocusRefRate(std::log(1-cleanLocusIndelRate));

    double logLhood(0.);
    for (const auto& obs : observations)
    {
        const double noisyMix(getObsLogLhood(logHomPrior, logHetPrior, logAltHetPrior, logNoIndelPrior,
                                             logInsertErrorRate, logDeleteErrorRate, logNoIndelRefRate, obs));
        const double cleanMix(getObsLogLhood(logHomPrior, logHetPrior, logAltHetPrior, logNoIndelPrior,
                                             logCleanLocusIndelRate, logCleanLocusIndelRate, logCleanLocusRefRate, obs));

        const double mix(log_sum(logCleanLocusRate+cleanMix,logNoisyLocusRate+noisyMix));

#ifdef DEBUG_MODEL3
        log_os << "MODEL3: loghood obs: noisy/clean/mix/delta: " << noisyMix << " " << cleanMix << " " << mix << " " << (mix*obs.repeatCount) << "\n";
#endif

        logLhood += (mix*obs.observationCount);
    }

#ifdef DEBUG_MODEL3
    log_os << "MODEL3: loghood output:" << logLhood << "\n";
#endif

    return logLhood;
}

struct error_minfunc_model3 : public codemin::minfunc_interface<double>
{
    explicit
    error_minfunc_model3(
        const std::vector<ExportedIndelObservations>& observations, const double theta,
        const bool isLockTheta = false)
        : defaultLogTheta(theta),
          _obs(observations),
          _isLockTheta(isLockTheta)
    {}
    unsigned dim() const override
    {
        return (_isLockTheta ? (MIN_PARAMS3::SIZE-1) : MIN_PARAMS3::SIZE);
    }

    double val(const double* in) override
    {
        argToParameters(in,_params);
        return -contextLogLhood(_obs,
                                _params[MIN_PARAMS3::LN_INSERT_ERROR_RATE],
                                _params[MIN_PARAMS3::LN_DELETE_ERROR_RATE],
                                _params[MIN_PARAMS3::LN_NOISY_LOCUS_RATE],
                                (_isLockTheta ? defaultLogTheta : _params[MIN_PARAMS3::LN_THETA]));
    }

    /// normalize the minimization values back to usable parameters
    ///
    /// most values are not valid on [-inf,inf] -- the minimizer doesn't
    /// know this. here is where we fill in the gap:
    ///
    static
    void
    argToParameters(
        const double* in,
        double* out)
    {
        auto rateSmoother = [](double a) -> double
        {
            static const double triggerVal(1e-3);
            static const double logTriggerVal(std::log(triggerVal));
            if (a>logTriggerVal)
            {
                a = std::log1p(a-logTriggerVal) + logTriggerVal;
            }
            return (a>maxLogRate ? maxLogRate-std::abs(a-maxLogRate) : a);
        };

        auto locusRateSmoother = [](double a) -> double
        {
            static const double triggerVal(0.8);
            static const double logTriggerVal(std::log(triggerVal));
            if (a>logTriggerVal)
            {
                a = std::log1p(a-logTriggerVal) + logTriggerVal;
            }
            return (a>maxLogLocusRate ? maxLogLocusRate-std::abs(a-maxLogLocusRate) : a);
        };

        // A lot of conditioning is required to keep the model from winding
        // theta around zero and getting confused, here we start applying a
        // second log to the delta above triggerTheta, and finally put a hard stop
        // at logMaxTheta -- hard stops are obviously bad b/c the model can get lost
        // on the flat plane even if the ML value is well below this limit, but
        // in practice this is such a ridiculously high value for theta, that
        // I don't see the model getting trapped.
        auto thetaSmoother = [](double a) -> double
        {
            static const double triggerVal(1e-3);
            static const double logTriggerVal(std::log(triggerVal));

            if (a>logTriggerVal)
            {
                a = std::log1p(a-logTriggerVal) + logTriggerVal;
            }
            return (a>maxLogTheta ? maxLogTheta-std::abs(a-maxLogTheta) : a);
        };

        for (unsigned paramIndex(MIN_PARAMS3::LN_INSERT_ERROR_RATE); paramIndex<MIN_PARAMS3::LN_NOISY_LOCUS_RATE; ++paramIndex)
        {
            out[paramIndex] = rateSmoother(in[paramIndex]);
        }
        out[MIN_PARAMS3::LN_NOISY_LOCUS_RATE] = locusRateSmoother(in[MIN_PARAMS3::LN_NOISY_LOCUS_RATE]);
        out[MIN_PARAMS3::LN_THETA] = thetaSmoother(in[MIN_PARAMS3::LN_THETA]);
    }

#if 0
    // this should help in theory, but in practice the minimizer is more likely to get stuck
    bool
    is_val_computable(
        const double* in) override
    {
        if (in[MIN_PARAMS3::LN_INSERT_ERROR_RATE]>maxLogRate) return false;
        if (in[MIN_PARAMS3::LN_DELETE_ERROR_RATE]>maxLogRate) return false;
        if (in[MIN_PARAMS3::LN_NOISY_LOCUS_RATE]>maxLogLocusRate) return false;
        if (in[MIN_PARAMS3::LN_THETA]>maxLogTheta) return false;
        return true;
    }
#endif

    const double defaultLogTheta;
    static const double initLogTheta;
    static const double maxLogTheta;
    static const double maxLogRate;
    static const double maxLogLocusRate;

private:
    const std::vector<ExportedIndelObservations>& _obs;
    bool _isLockTheta;
    double _params[MIN_PARAMS3::SIZE];
};

const double error_minfunc_model3::initLogTheta = std::log(1e-04);
const double error_minfunc_model3::maxLogTheta = std::log(0.4);
const double error_minfunc_model3::maxLogRate = std::log(0.5);
const double error_minfunc_model3::maxLogLocusRate = std::log(1.0);

static
bool
computeExtendedContext(
    const bool isLockTheta,
    const double logTheta,
    const IndelErrorData& data,
    double normalizedParams[MIN_PARAMS3::SIZE])
{
    bool paramsAcceptable = true;
    std::vector<ExportedIndelObservations> observations;
    data.exportObservations(observations);
    // initialize conjugate direction minimizer settings and minimize lhood...
    //
    double minParams[MIN_PARAMS3::SIZE];

    {
        unsigned iter;
        double x_all_loghood;
        static const double line_tol(1e-10);
        static const double end_tol(1e-10);
        static const unsigned max_iter(40);

        error_minfunc_model3 errFunc(observations, logTheta, isLockTheta);
        // initialize parameter search
        minParams[MIN_PARAMS3::LN_INSERT_ERROR_RATE] = std::log(1e-3);
        minParams[MIN_PARAMS3::LN_DELETE_ERROR_RATE] = std::log(1e-3);
        minParams[MIN_PARAMS3::LN_NOISY_LOCUS_RATE] = std::log(0.4);
        minParams[MIN_PARAMS3::LN_THETA] = errFunc.defaultLogTheta;

        static const unsigned SIZE2(MIN_PARAMS3::SIZE*MIN_PARAMS3::SIZE);
        double conjDir[SIZE2];

        std::fill(conjDir,conjDir+SIZE2,0.);
        const unsigned dim(isLockTheta ? MIN_PARAMS3::SIZE-1 : MIN_PARAMS3::SIZE);
        for (unsigned i(0); i<dim; ++i)
        {
            conjDir[i*(dim+1)] = 0.0005;
        }

        double start_tol(end_tol);
        double final_dlh;


        codemin::minimize_conj_direction(minParams,conjDir,errFunc,start_tol,end_tol,line_tol,
                                         x_all_loghood,iter,final_dlh,max_iter);

        if (max_iter == iter)
        {
            paramsAcceptable = false;
        }
    }

    error_minfunc_model3::argToParameters(minParams,normalizedParams);

    return paramsAcceptable;
}



static
AdaptiveIndelErrorModelLogParams
estimateModelParams(
    const SequenceErrorCounts& counts,
    const IndelErrorContext context,
    const double logTheta)
{
    // setup the optimizer settings to the model assumption
    const bool isLockTheta = true;
    double normalizedParams[MIN_PARAMS3::SIZE];

    AdaptiveIndelErrorModelLogParams estimatedParams;
    auto contextIt = counts.getIndelCounts().find(context);
    if (contextIt != counts.getIndelCounts().end())
    {

        const auto& data(contextIt->second);

        estimatedParams.paramsAcceptable = computeExtendedContext(isLockTheta, logTheta, data, normalizedParams);

        estimatedParams.logErrorRate = (normalizedParams[MIN_PARAMS3::LN_INSERT_ERROR_RATE] +
                                        normalizedParams[MIN_PARAMS3::LN_DELETE_ERROR_RATE]) / 2;
        estimatedParams.logNoisyLocusRate = normalizedParams[MIN_PARAMS3::LN_NOISY_LOCUS_RATE];
    }
    return estimatedParams;
}

IndelModelProduction::
IndelModelProduction(
    const SequenceErrorCounts& counts,
    const std::string& thetaFilename,
    const std::string& outputFilename):
    _counts(counts),
    _outputFilename(outputFilename)
{
    if (thetaFilename.empty())
    {
        using namespace illumina::common;

        std::ostringstream oss;
        oss << "Theta file name cannot be empty\n";
        BOOST_THROW_EXCEPTION(LogicException(oss.str()));
    }

    _thetas = IndelErrorModelJson::deserializeTheta(thetaFilename);

}

void
IndelModelProduction::
estimateIndelErrorRates()
{
    const auto lowRepeatCount = AdaptiveIndelErrorModel::lowRepeatCount;
    assert(_repeatPatterns.size() == _maxRepeatCounts.size());

    for (unsigned repeatPatternIndex = 0; repeatPatternIndex < _repeatPatterns.size(); repeatPatternIndex++)
    {
        auto repeatPatternSize = _repeatPatterns[repeatPatternIndex];
        auto theta = _thetas[repeatPatternSize];
        assert(theta.size() >= *std::max_element(_maxRepeatCounts.begin(), _maxRepeatCounts.end()));

        // estimate low repeat count params
        const auto highRepeatCount = _maxRepeatCounts[repeatPatternIndex];
        IndelErrorContext lowCountContext(repeatPatternSize, lowRepeatCount);
        log_os << "INFO: computing rates for context: " << lowCountContext << "\n";
        const auto lowLogParams = estimateModelParams(_counts, lowCountContext, std::log(theta[lowRepeatCount - 1]));
        // estimate high repeat count params
        IndelErrorContext highCountContext(repeatPatternSize, highRepeatCount);
        log_os << "INFO: computing rates for context: " << highCountContext << "\n";
        const auto highLogParams = estimateModelParams(_counts, highCountContext, std::log(theta[highRepeatCount - 1]));

        _adaptiveIndelErrorModels.push_back(AdaptiveIndelErrorModel(repeatPatternSize,
                                                                    highRepeatCount,
                                                                    lowLogParams,
                                                                    highLogParams));
        if (!lowLogParams.paramsAcceptable || !highLogParams.paramsAcceptable)
        {
            _isEstimationAcceptable = false;
        }
    }

    // estimate error rate for the non-STR context
    const unsigned nonSTRRepeatPatternSize(1);
    const unsigned nonSTRRepeatCount(1);
    IndelErrorContext targetContext(nonSTRRepeatPatternSize, nonSTRRepeatCount);
    const auto nonSTRTheta = _thetas.at(nonSTRRepeatPatternSize)[0];
    log_os << "INFO: computing rates for context: " << targetContext << "\n";
    _nonSTRModelParams = estimateModelParams(_counts, targetContext, std::log(nonSTRTheta));
    _isEstimated = true;
}

void
IndelModelProduction::exportModel() const
{
    IndelErrorModelJson indelModelJson(_counts.getSampleName());
    // add the non-STR params to all contexts with repeat count 1
    // this will show up as valid contexts during variant calling so we need to fill in these gaps
    for (unsigned repeatPatternIndex = 0; repeatPatternIndex < _repeatPatterns.size(); repeatPatternIndex++)
    {
        indelModelJson.addMotif(_repeatPatterns[repeatPatternIndex], 1, std::exp(_nonSTRModelParams.logErrorRate), std::exp(_nonSTRModelParams.logNoisyLocusRate));
    }

    // add motif to json for all contexts
    for (unsigned repeatPatternIndex = 0; repeatPatternIndex < _repeatPatterns.size(); repeatPatternIndex++)
    {
        auto errorModel = _adaptiveIndelErrorModels[repeatPatternIndex];

        for (unsigned repeatCount = errorModel.lowRepeatCount; repeatCount <=errorModel.highRepeatCount(); repeatCount++)
        {
            indelModelJson.addMotif(errorModel.repeatPatternSize(),
                                    repeatCount,
                                    errorModel.errorRate(repeatCount),
                                    errorModel.noisyLocusRate(repeatCount));
        }
    }

    static const bool isStatic(false);
    indelModelJson.serializeIndelErrorModel(indelModelJson.getSampleName(),
                                            indelModelJson.generateMotifsNode(),
                                            isStatic,
                                            _outputFilename);
}

void
IndelModelProduction::exportModelUsingInputJson(const std::string& jsonFilename) const
{
    std::string jsonString;
    Json::Value root;
    {
        std::ifstream ifs(jsonFilename, std::ifstream::binary);
        std::stringstream buffer;
        buffer << ifs.rdbuf();
        jsonString = buffer.str();
    }
    Json::Reader reader;
    if (!reader.parse(jsonString, root))
    {
        using namespace illumina::common;

        std::ostringstream oss;
        oss << "Failed to parse JSON " << jsonFilename << " " << reader.getFormattedErrorMessages() << "'\n";
        BOOST_THROW_EXCEPTION(LogicException(oss.str()));
    }

    Json::Value samples = root["sample"];
    if (samples.isNull() || samples.empty())
    {
        using namespace illumina::common;
        std::ostringstream oss;
        oss << "ERROR: no samples in indel error model file '" << jsonFilename << "'\n";
        BOOST_THROW_EXCEPTION(LogicException(oss.str()));
    }
    else if (samples.size() > 1)
    {
        using namespace illumina::common;
        std::ostringstream oss;
        oss << "ERROR: multiple samples in indel error model file '" << jsonFilename << "'\n";
        BOOST_THROW_EXCEPTION(LogicException(oss.str()));
    }

    Json::Value motifs = samples[0]["motif"];
    if (motifs.isNull() || motifs.empty())
    {
        using namespace illumina::common;
        std::ostringstream oss;
        oss << "ERROR: no motifs in indel error model file '" << jsonFilename << "'\n";
        BOOST_THROW_EXCEPTION(LogicException(oss.str()));
    }

    static const bool isStatic(true);
    IndelErrorModelJson::serializeIndelErrorModel(_counts.getSampleName(), motifs, isStatic, _outputFilename);
}

bool
IndelModelProduction::
checkEstimatedModel() const
{
    if (!_isEstimated)
    {
        using namespace illumina::common;
        std::ostringstream oss;
        oss << "ERROR: indel error model has not been estimated '\n";
        BOOST_THROW_EXCEPTION(LogicException(oss.str()));
    }

    // check estimated STR params
    for (const auto& model : _adaptiveIndelErrorModels)
    {
        if (!isValidErrorRate(model.errorRate(model.lowRepeatCount)) ||
            !isValidErrorRate(model.errorRate(model.highRepeatCount())))
        {
            return false;
        }

    }

    // check non-STR params and the _isEstimationAcceptable flag
    if (!isValidErrorRate(std::exp(_nonSTRModelParams.logErrorRate)) || !_isEstimationAcceptable)
    {
        return false;
    }

    return true;
}

bool IndelModelProduction::
isValidErrorRate(
    const double errorRate) const
{
    return !(std::isnan(errorRate) ||
             std::isinf(errorRate) ||
             (errorRate >= _maxErrorRate) ||
             (errorRate <= _minErrorRate));
}

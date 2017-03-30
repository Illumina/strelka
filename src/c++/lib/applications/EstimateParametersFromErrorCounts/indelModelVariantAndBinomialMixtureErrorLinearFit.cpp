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

///
/// \author Chris Saunders
///

#include "blt_util/log.hh"
#include "blt_util/math_util.hh"
#include "blt_util/prob_util.hh"
#include "indelModelVariantAndBinomialMixtureErrorLinearFit.hh"

//#define CODEMIN_DEBUG
#define CODEMIN_USE_BOOST
#include "minimize_conj_direction.h"

#include <iomanip>
#include <iostream>

namespace
{

namespace MIN_PARAMS3_SIMPLE
{
enum index_t
{
    LN_INDEL_ERROR_RATE_MIN,
    LN_INDEL_ERROR_RATE_MAX,
    LN_NOISY_LOCUS_RATE_MIN,
    LN_NOISY_LOCUS_RATE_MAX,
    LN_THETA,
    SIZE
};
}


//#define DEBUG_MODEL3



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
        const double logInsertErrorRateSlope,
        const double logInsertErrorRateIntercept,
        const double logDeleteErrorRateSlope,
        const double logDeleteErrorRateIntercept,
        const double logNoisyLocusRateSlope,
        const double logNoisyLocusRateIntercept,
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


    double logLikelihood(0.);
    for (const auto& obs : observations)
    {
        const auto repeatCount = obs.repeatCount;

        const auto logInsertErrorRate = repeatCount * logInsertErrorRateSlope + logInsertErrorRateIntercept;
        const auto logDeleteErrorRate = repeatCount * logDeleteErrorRateSlope + logDeleteErrorRateIntercept;
        const auto logNoisyLocusRate = repeatCount * logNoisyLocusRateSlope + logNoisyLocusRateIntercept;

        const double logNoIndelRefRate(std::log(1-std::exp(logInsertErrorRate))+std::log(1-std::exp(logDeleteErrorRate)));

        const double logCleanLocusRate(std::log(1-std::exp(logNoisyLocusRate)));

        static const double cleanLocusIndelRate(1e-8);
        static const double logCleanLocusIndelRate(std::log(cleanLocusIndelRate));
        static const double logCleanLocusRefRate(std::log(1-cleanLocusIndelRate));



        const double noisyMix(getObsLogLhood(logHomPrior, logHetPrior, logAltHetPrior, logNoIndelPrior,
                                             logInsertErrorRate, logDeleteErrorRate, logNoIndelRefRate, obs));
        const double cleanMix(getObsLogLhood(logHomPrior, logHetPrior, logAltHetPrior, logNoIndelPrior,
                                             logCleanLocusIndelRate, logCleanLocusIndelRate, logCleanLocusRefRate, obs));

        const double mix(log_sum(logCleanLocusRate+cleanMix,logNoisyLocusRate+noisyMix));

#ifdef DEBUG_MODEL3
        log_os << "MODEL3: loghood obs: noisy/clean/mix/delta: " << noisyMix << " " << cleanMix << " " << mix << " " << (mix*obs.repeatCount) << "\n";
#endif

        logLikelihood += (mix*obs.observationCount);
    }

#ifdef DEBUG_MODEL3
    log_os << "MODEL3: loghood output:" << logLhood << "\n";
#endif

    return logLikelihood;
}


struct errorMinfuncModel3Simple : public codemin::minfunc_interface<double>
{
    explicit
    errorMinfuncModel3Simple(
            const std::vector<ExportedIndelObservations>& observations,
            const unsigned minSTRLength,
            const unsigned maxSTRLength,
            const bool isLockTheta = true)
            : _obs(observations), _minSTRLength(minSTRLength), _maxSTRLength(maxSTRLength), _isLockTheta(isLockTheta)
    {}

    unsigned dim() const override
    {
        return (MIN_PARAMS3_SIMPLE::SIZE-1);
    }

    double val(const double* in) override
    {
        argToParameters(in,_params);

        const auto diffX = (double)(_maxSTRLength - _minSTRLength);
        const double logInsertErrorRateSlope = (_params[MIN_PARAMS3_SIMPLE::LN_INDEL_ERROR_RATE_MAX] - _params[MIN_PARAMS3_SIMPLE::LN_INDEL_ERROR_RATE_MIN])/diffX;
        const double logInsertErrorRateIntercept = _params[MIN_PARAMS3_SIMPLE::LN_INDEL_ERROR_RATE_MIN] - logInsertErrorRateSlope * _minSTRLength;
        const double logDeleteErrorRateSlope = (_params[MIN_PARAMS3_SIMPLE::LN_INDEL_ERROR_RATE_MAX] - _params[MIN_PARAMS3_SIMPLE::LN_INDEL_ERROR_RATE_MIN])/diffX;
        const double logDeleteErrorRateIntercept = _params[MIN_PARAMS3_SIMPLE::LN_INDEL_ERROR_RATE_MIN] - logDeleteErrorRateSlope * _minSTRLength;
        const double logNoisyLocusRateSlope = (_params[MIN_PARAMS3_SIMPLE::LN_NOISY_LOCUS_RATE_MAX] - _params[MIN_PARAMS3_SIMPLE::LN_NOISY_LOCUS_RATE_MIN])/diffX;
        const double logNoisyLocusRateIntercept = _params[MIN_PARAMS3_SIMPLE::LN_NOISY_LOCUS_RATE_MIN] - logNoisyLocusRateSlope * _minSTRLength;

        return -contextLogLhood(_obs,
                                logInsertErrorRateSlope,
                                logInsertErrorRateIntercept,
                                logDeleteErrorRateSlope,
                                logDeleteErrorRateIntercept,
                                logNoisyLocusRateSlope,
                                logNoisyLocusRateIntercept,
                                defaultLogTheta);
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
                a = std::log(1+(a-logTriggerVal)) + logTriggerVal;
            }
            return (a>maxLogRate ? maxLogRate-std::abs(a-maxLogRate) : a);
        };

        auto locusRateSmoother = [](double a) -> double
        {
            static const double triggerVal(0.8);
            static const double logTriggerVal(std::log(triggerVal));
            if (a>logTriggerVal)
            {
                a = std::log(1+(a-logTriggerVal)) + logTriggerVal;
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
                a = std::log(1+(a-logTriggerVal)) + logTriggerVal;
            }
            return (a>maxLogTheta ? maxLogTheta-std::abs(a-maxLogTheta) : a);
        };

        // this shouldn't really work yet
        for (unsigned paramIndex(MIN_PARAMS3_SIMPLE::LN_INDEL_ERROR_RATE_MIN); paramIndex<=MIN_PARAMS3_SIMPLE::LN_INDEL_ERROR_RATE_MAX; ++paramIndex)
        {
            out[paramIndex] = rateSmoother(in[paramIndex]);
        }
        for (unsigned paramIndex(MIN_PARAMS3_SIMPLE::LN_NOISY_LOCUS_RATE_MIN); paramIndex<=MIN_PARAMS3_SIMPLE::LN_NOISY_LOCUS_RATE_MAX; ++paramIndex)
        {
            out[paramIndex] = locusRateSmoother(in[paramIndex]);
        }

        out[MIN_PARAMS3_SIMPLE::LN_THETA] = thetaSmoother(in[MIN_PARAMS3_SIMPLE::LN_THETA]);
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

    static const double defaultLogTheta;
    static const double maxLogTheta;
    static const double maxLogRate;
    static const double maxLogLocusRate;

private:
    const std::vector<ExportedIndelObservations>& _obs;
    const unsigned _minSTRLength;
    const unsigned _maxSTRLength;
    bool _isLockTheta;
    double _params[MIN_PARAMS3_SIMPLE::SIZE];
};

const double errorMinfuncModel3Simple::defaultLogTheta = std::log(1e-4);
const double errorMinfuncModel3Simple::maxLogTheta = std::log(0.4);
const double errorMinfuncModel3Simple::maxLogRate = std::log(0.5);
const double errorMinfuncModel3Simple::maxLogLocusRate = std::log(1.0);



struct SignalGroupTotal
{
    double ref = 0;
    double alt = 0;
    double locus = 0;
};



static
void
getAltSigTotal(
        const std::vector<ExportedIndelObservations>& observations,
        const unsigned altBeginIndex,
        const unsigned altEndIndex,
        SignalGroupTotal& sigTotal)
{
    for (const ExportedIndelObservations& obs : observations)
    {
        unsigned totalAltObservations(0);
        for (unsigned altIndex(altBeginIndex); altIndex<altEndIndex; ++altIndex)
        {
            totalAltObservations += obs.altObservations[altIndex];
        }

        sigTotal.ref += (obs.refObservations*obs.observationCount);
        sigTotal.alt += (totalAltObservations*obs.observationCount);
        sigTotal.locus += obs.observationCount;
    }
}



static
void
reportIndelErrorRateSet(
        const IndelErrorContext& context,
        const char* extendedContextTag,
        unsigned iter,
        const double loghood,
        const double indelErrorRate,
        const double theta,
        const double noisyLocusRate,
        std::ostream& os)
{
    static const std::string sep(", ");

    os << std::setprecision(10);
    os << context << "_" << extendedContextTag
       << sep << iter
       << sep << loghood
       << sep << indelErrorRate
       << sep << theta
       << sep << noisyLocusRate
       << "\n";
}



static
void
estimateParameters(
        const std::vector<ExportedIndelObservations>& observations,
        unsigned minSTRLength,
        unsigned maxSTRLength,
        std::ostream& os)
{
    // Get summary counts for QC purposes. Note these are unrelated to minimization or model:
    SignalGroupTotal sigInsertTotal;
    getAltSigTotal(observations, INDEL_SIGNAL_TYPE::INSERT_1, INDEL_SIGNAL_TYPE::DELETE_1, sigInsertTotal);

    SignalGroupTotal sigDeleteTotal;
    getAltSigTotal(observations, INDEL_SIGNAL_TYPE::DELETE_1, INDEL_SIGNAL_TYPE::SIZE, sigDeleteTotal);


    // initialize conjugate direction minimizer settings and minimize lhood...
    //
    double minParams[MIN_PARAMS3_SIMPLE::SIZE];

    unsigned iter;
    double x_all_loghood;
    {
        static const double line_tol(1e-10);
        static const double end_tol(1e-10);
        static const unsigned max_iter(40);

        // initialize parameter search
        minParams[MIN_PARAMS3_SIMPLE::LN_INDEL_ERROR_RATE_MIN] = std::log(1e-3);
        minParams[MIN_PARAMS3_SIMPLE::LN_INDEL_ERROR_RATE_MAX] = std::log(1e-3);
        minParams[MIN_PARAMS3_SIMPLE::LN_NOISY_LOCUS_RATE_MIN] = std::log(1e-3);
        minParams[MIN_PARAMS3_SIMPLE::LN_NOISY_LOCUS_RATE_MAX] = std::log(1e-3);
        minParams[MIN_PARAMS3_SIMPLE::LN_THETA] = errorMinfuncModel3Simple::defaultLogTheta;

        static const unsigned SIZE2(MIN_PARAMS3_SIMPLE::SIZE*MIN_PARAMS3_SIMPLE::SIZE);
        double conjDir[SIZE2];

        std::fill(conjDir,conjDir+SIZE2,0.);
        const unsigned dim(MIN_PARAMS3_SIMPLE::SIZE-1);
        for (unsigned i(0); i<dim; ++i)
        {
            conjDir[i*(dim+1)] = 0.0005;
        }

        double start_tol(end_tol);
        double final_dlh;
        errorMinfuncModel3Simple errFunc(observations, minSTRLength, maxSTRLength, true);

        codemin::minimize_conj_direction(minParams,conjDir,errFunc,start_tol,end_tol,line_tol,
                                         x_all_loghood,iter,final_dlh,max_iter);
    }

    // We will not need this at the end... this is to conveniently compare to the previous results
    {

        double normalizedParams[MIN_PARAMS3_SIMPLE::SIZE];
        errorMinfuncModel3Simple::argToParameters(minParams,normalizedParams);
        const double theta(std::exp(normalizedParams[MIN_PARAMS3_SIMPLE::LN_THETA]));
        const auto diffX = (maxSTRLength - minSTRLength);
        const double logInsertErrorRateSlope = (normalizedParams[MIN_PARAMS3_SIMPLE::LN_INDEL_ERROR_RATE_MAX] - normalizedParams[MIN_PARAMS3_SIMPLE::LN_INDEL_ERROR_RATE_MIN])/diffX;
        const double logInsertErrorRateIntercept = normalizedParams[MIN_PARAMS3_SIMPLE::LN_INDEL_ERROR_RATE_MIN] - logInsertErrorRateSlope * minSTRLength;
        const double logDeleteErrorRateSlope = (normalizedParams[MIN_PARAMS3_SIMPLE::LN_INDEL_ERROR_RATE_MAX] - normalizedParams[MIN_PARAMS3_SIMPLE::LN_INDEL_ERROR_RATE_MIN])/diffX;
        const double logDeleteErrorRateIntercept = normalizedParams[MIN_PARAMS3_SIMPLE::LN_INDEL_ERROR_RATE_MIN] - logDeleteErrorRateSlope * minSTRLength;
        const double logNoisyLocusRateSlope = (normalizedParams[MIN_PARAMS3_SIMPLE::LN_NOISY_LOCUS_RATE_MAX] - normalizedParams[MIN_PARAMS3_SIMPLE::LN_NOISY_LOCUS_RATE_MIN])/diffX;
        const double logNoisyLocusRateIntercept = normalizedParams[MIN_PARAMS3_SIMPLE::LN_NOISY_LOCUS_RATE_MIN] - logNoisyLocusRateSlope * minSTRLength;

        log_os << "INFO: logInsertErrorRateSlope: " << logInsertErrorRateSlope << "\n";
        log_os << "INFO: logInsertErrorRateIntercept: " << logInsertErrorRateIntercept << "\n";

        log_os << "INFO: logDeleteErrorRateSlope: " << logDeleteErrorRateSlope << "\n";
        log_os << "INFO: logDeleteErrorRateIntercept: " << logDeleteErrorRateIntercept << "\n";

        log_os << "INFO: logNoisyLocusRateSlope: " << logNoisyLocusRateSlope << "\n";
        log_os << "INFO: logNoisyLocusRateIntercept: " << logNoisyLocusRateIntercept << "\n";

        for (unsigned repeatCount = minSTRLength; repeatCount <= maxSTRLength; repeatCount++)
        {

            const IndelErrorContext context(observations.front().repeatPatternSize, repeatCount);

            const auto insertErrorRate(std::exp(repeatCount * logInsertErrorRateSlope + logInsertErrorRateIntercept));
            const auto deleteErrorRate(std::exp(repeatCount * logDeleteErrorRateSlope + logDeleteErrorRateIntercept));
            const auto noisyLocusRate(std::exp(repeatCount * logNoisyLocusRateSlope + logNoisyLocusRateIntercept));



            reportIndelErrorRateSet(context, "I", iter, -x_all_loghood, insertErrorRate, theta,
                                    noisyLocusRate, os);
            reportIndelErrorRateSet(context, "D", iter, -x_all_loghood, deleteErrorRate, theta,
                                    noisyLocusRate, os);
        }
    }
}

}



void
indelModelVariantAndBinomialMixtureErrorLinearFit(
        const SequenceErrorCounts &counts)
{
    const std::vector<unsigned> referenceSTRPatternSizeVector = {1,2};
    std::ostream& ros(std::cout);

    const std::vector<unsigned> minSTRLengthVector = {2, 2};
    const std::vector<unsigned> maxSTRLengthVector = {16, 8};

    assert(referenceSTRPatternSizeVector.size() == minSTRLengthVector.size() && referenceSTRPatternSizeVector.size() == maxSTRLengthVector.size());

    ros << "context, iter, lhood, errorRate, theta, noisyLocusRate\n";

    // only consider repeat counts > 1 for now
    for (unsigned indexSTR = 0; indexSTR < referenceSTRPatternSizeVector.size(); indexSTR++)
    {
        const auto referenceSTRPatternSize = referenceSTRPatternSizeVector[indexSTR];
        const auto minSTRLength = minSTRLengthVector[indexSTR];
        const auto maxSTRLength = maxSTRLengthVector[indexSTR];

        assert(maxSTRLength > minSTRLength); // has to be strictly larger

        std::vector<ExportedIndelObservations> observations;
        for (const auto &contextInfo : counts.getIndelCounts())
        {
            if (contextInfo.first.getRepeatPatternSize() == referenceSTRPatternSize &&
                    contextInfo.first.getRepeatCount() >= minSTRLength &&
                    contextInfo.first.getRepeatCount() <= maxSTRLength)
            {
                const auto& context(contextInfo.first);
                const auto& data(contextInfo.second);
                data.addToObservations(context, observations);
            }
        }

        if (observations.empty())
        {
            log_os << "ERROR: No errors observed for pattern size: " << referenceSTRPatternSize << "\n";
            return;
        }
        log_os << "INFO: computing rates for all contexts with repeatPatternSize == " << referenceSTRPatternSize << " and repeatCount " << minSTRLength << " " << maxSTRLength << "\n";
        estimateParameters(observations, minSTRLength, maxSTRLength, ros);
    }
}

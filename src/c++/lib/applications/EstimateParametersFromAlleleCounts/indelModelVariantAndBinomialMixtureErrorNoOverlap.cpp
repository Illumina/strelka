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

/// \file
/// \author Chris Saunders
///

#include "indelModelVariantAndBinomialMixtureErrorNoOverlap.hh"

#include "blt_util/log.hh"
#include "blt_util/logSumUtil.hh"
#include "blt_util/prob_util.hh"

//#define CODEMIN_DEBUG
#define CODEMIN_USE_BOOST
#include "minimize_conj_direction.h"

#include <iomanip>
#include <iostream>

namespace
{

namespace MIN_PARAMS3
{
enum index_t
{
    LN_INDEL_ERROR_RATE,
    LN_NOISY_LOCUS_RATE,
    LN_THETA,
    SIZE
};
}


//#define DEBUG_MODEL3


using namespace IndelCounts;


static const double cleanLocusIndelErrorRate(1e-8);



static
double
getObsLogLhood(
    const double logHomPrior,
    const double logHetPrior,
    const double logNoIndelPrior,
    const double logIndelErrorRate,
    const double logNoIndelRefRate,
    const bool isInsert,
    const IndelCounts::SingleSampleContextObservationInfoExportFormat& contextObservationInfo)
{
    static const double homAltRate(0.99);
    static const double hetAltRate(0.5);

    static const double logHomAltRate(std::log(homAltRate));
    static const double logHomRefRate(std::log(1.-homAltRate));
    static const double logHetRate(std::log(hetAltRate));

    unsigned totalIndelObservations(0);
    if (isInsert)
    {
        for (unsigned altIndex(INDEL_SIGNAL_TYPE::INSERT_1); altIndex<INDEL_SIGNAL_TYPE::DELETE_1; ++altIndex)
        {
            totalIndelObservations += contextObservationInfo.altObservations[altIndex];
        }
    }
    else
    {
        for (unsigned altIndex(INDEL_SIGNAL_TYPE::DELETE_1); altIndex<INDEL_SIGNAL_TYPE::SIZE; ++altIndex)
        {
            totalIndelObservations += contextObservationInfo.altObservations[altIndex];
        }
    }

    // get lhood of homref GT:
    const double noindel(logIndelErrorRate*totalIndelObservations +
                         logNoIndelRefRate*contextObservationInfo.refObservations);

    // get lhood of het and hom GT:
    const double het(logHetRate*(contextObservationInfo.refObservations+totalIndelObservations));

    const double hom(logHomAltRate*totalIndelObservations + logHomRefRate*contextObservationInfo.refObservations);

    return getLogSum(logHomPrior+hom, logHetPrior+het, logNoIndelPrior+noindel);
}



static
double
contextLogLhood(
    const SingleSampleContextDataExportFormat& exportedContextData,
    const double logIndelErrorRate,
    const double logNoisyLocusRate,
    const bool isInsert,
    const double logTheta)
{
#ifdef DEBUG_MODEL3
    log_os << "MODEL3: loghood input:"
           << " insert: " << std::exp(logIndelErrorRate)
           << " noise: " << std::exp(logNoisyLocusRate)
           << " theta: " << std::exp(logTheta)
           << "\n";
#endif

    static const double log2(std::log(2));
    const double logHomPrior(logTheta-log2);
    const double logHetPrior(logTheta);
    const double theta(std::exp(logTheta));
    const double logNoIndelPrior(std::log(1-(theta*3./2.)));

    const double logNoIndelRefRate(std::log1p(-std::exp(logIndelErrorRate)));

    const double logCleanLocusRate(std::log(1-std::exp(logNoisyLocusRate)));

    static const double logCleanLocusIndelRate(std::log(cleanLocusIndelErrorRate));
    static const double logCleanLocusRefRate(std::log(1-cleanLocusIndelErrorRate));

    double logLhood(0.);
    for (const auto& contextObservationInfo : exportedContextData.data)
    {
        const double noisyMix(getObsLogLhood(logHomPrior, logHetPrior, logNoIndelPrior,
                                             logIndelErrorRate, logNoIndelRefRate, isInsert, contextObservationInfo));
        const double cleanMix(getObsLogLhood(logHomPrior, logHetPrior, logNoIndelPrior,
                                             logCleanLocusIndelRate, logCleanLocusRefRate, isInsert, contextObservationInfo));

        const double mix(getLogSum(logCleanLocusRate+cleanMix, logNoisyLocusRate+noisyMix));

#ifdef DEBUG_MODEL3
        log_os << "MODEL3: loghood obs: noisy/clean/mix/delta: " << noisyMix << " " << cleanMix << " " << mix << " " << (mix*obs.repeatCount) << "\n";
#endif

        logLhood += (mix*contextObservationInfo.contextInstanceCount);
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
        const SingleSampleContextDataExportFormat& exportedContextData,
        const bool isInsert,
        const bool isLockTheta = false)
        : _exportedContextData(exportedContextData),
          _isInsert(isInsert),
          _isLockTheta(isLockTheta)
    {}

    unsigned dim() const override
    {
        return (_isLockTheta ? (MIN_PARAMS3::SIZE-1) : MIN_PARAMS3::SIZE);
    }

    double val(const double* in) override
    {
        argToParameters(in,_params);
        return -contextLogLhood(_exportedContextData,
                                _params[MIN_PARAMS3::LN_INDEL_ERROR_RATE],
                                _params[MIN_PARAMS3::LN_NOISY_LOCUS_RATE],
                                _isInsert,
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

        for (unsigned paramIndex(MIN_PARAMS3::LN_INDEL_ERROR_RATE); paramIndex<MIN_PARAMS3::LN_NOISY_LOCUS_RATE; ++paramIndex)
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

    static const double defaultLogTheta;
    static const double maxLogTheta;
    static const double maxLogRate;
    static const double maxLogLocusRate;

private:
    const SingleSampleContextDataExportFormat& _exportedContextData;
    bool _isInsert;
    bool _isLockTheta;
    double _params[MIN_PARAMS3::SIZE];
};

const double error_minfunc_model3::defaultLogTheta = std::log(1e-4);
const double error_minfunc_model3::maxLogTheta = std::log(0.4);
const double error_minfunc_model3::maxLogRate = std::log(0.5);
const double error_minfunc_model3::maxLogLocusRate = std::log(1.0);



struct SignalGroupTotal
{
    double ref = 0;
    double alt = 0;
    double locus = 0;
};



static
void
getAltSigTotal(
    const SingleSampleContextDataExportFormat& exportedContextData,
    const unsigned altBeginIndex,
    const unsigned altEndIndex,
    SignalGroupTotal& sigTotal)
{
    for (const auto& contextObservationInfo : exportedContextData.data)
    {
        unsigned totalAltObservations(0);
        for (unsigned altIndex(altBeginIndex); altIndex<altEndIndex; ++altIndex)
        {
            totalAltObservations += contextObservationInfo.altObservations[altIndex];
        }

        sigTotal.ref += (contextObservationInfo.refObservations*contextObservationInfo.contextInstanceCount);
        sigTotal.alt += (totalAltObservations*contextObservationInfo.contextInstanceCount);
        sigTotal.locus += contextObservationInfo.contextInstanceCount;
    }
}



static
void
reportIndelErrorRateSet(
    const Context& context,
    const char* extendedContextTag,
    const SignalGroupTotal& sigTotal,
    const ContextData& contextData,
    unsigned iter,
    const double loghood,
    const double noisyLocusIndelErrorRate,
    const double noisyLocusRate,
    const double theta,
    std::ostream& os)
{
    static const std::string sep(", ");

    const double simpleIndelErrorRate(noisyLocusIndelErrorRate*noisyLocusRate + cleanLocusIndelErrorRate*(1-noisyLocusRate));

    os << std::setprecision(10);
    os << context << "_" << extendedContextTag
       << sep << contextData.excludedRegionSkipped
       << sep << (sigTotal.locus + contextData.depthSkipped)
       << sep << sigTotal.locus
       << sep << sigTotal.ref
       << sep << sigTotal.alt
       << sep << iter
       << sep << loghood
       << sep << noisyLocusIndelErrorRate
       << sep << cleanLocusIndelErrorRate
       << sep << noisyLocusRate
       << sep << simpleIndelErrorRate
       << sep << theta
       << "\n";
}



static
void
reportExtendedContext(
    const bool isLockTheta,
    const Context& context,
    const SingleSampleContextDataExportFormat& exportedContextData,
    const ContextData& contextData,
    std::ostream& os)
{
    // Get summary counts for QC purposes. Note these are unrelated to minimization or model:
    SignalGroupTotal sigInsertTotal;
    getAltSigTotal(exportedContextData, INDEL_SIGNAL_TYPE::INSERT_1, INDEL_SIGNAL_TYPE::DELETE_1, sigInsertTotal);

    SignalGroupTotal sigDeleteTotal;
    getAltSigTotal(exportedContextData, INDEL_SIGNAL_TYPE::DELETE_1, INDEL_SIGNAL_TYPE::SIZE, sigDeleteTotal);


    // initialize conjugate direction minimizer settings and minimize lhood...
    //
    for (unsigned indelTypeIndex(0); indelTypeIndex<2; ++indelTypeIndex)
    {
        const bool isInsert(indelTypeIndex==0);
        double minParams[MIN_PARAMS3::SIZE];

        unsigned iter;
        double x_all_loghood;
        {
            static const double line_tol(1e-10);
            static const double end_tol(1e-10);
            static const unsigned max_iter(40);

            // initialize parameter search
            minParams[MIN_PARAMS3::LN_INDEL_ERROR_RATE] = std::log(1e-3);
            minParams[MIN_PARAMS3::LN_NOISY_LOCUS_RATE] = std::log(0.4);
            minParams[MIN_PARAMS3::LN_THETA] = error_minfunc_model3::defaultLogTheta;

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
            error_minfunc_model3 errFunc(exportedContextData, isInsert, isLockTheta);

            codemin::minimize_conj_direction(minParams,conjDir,errFunc,start_tol,end_tol,line_tol,
                                             x_all_loghood,iter,final_dlh,max_iter);
        }

        // report:
        {
            double normalizedParams[MIN_PARAMS3::SIZE];
            error_minfunc_model3::argToParameters(minParams,normalizedParams);

            const double theta(std::exp(normalizedParams[MIN_PARAMS3::LN_THETA]));

            const double noisyLocusIndelErrorRate(std::exp(normalizedParams[MIN_PARAMS3::LN_INDEL_ERROR_RATE]));
            const double noisyLocusRate(std::exp(normalizedParams[MIN_PARAMS3::LN_NOISY_LOCUS_RATE]));
            const std::string tag(isInsert ? "I" : "D");
            reportIndelErrorRateSet(context, tag.c_str(), sigInsertTotal, contextData, iter, -x_all_loghood, noisyLocusIndelErrorRate, noisyLocusRate, theta, os);
        }
    }
}

}



void
indelModelVariantAndBinomialMixtureErrorNoOverlap(
    const SequenceAlleleCounts& counts)
{
    const bool isLockTheta(false);

    std::ostream& ros(std::cout);

    ros << "context, excludedLoci, nonExcludedLoci, usedLoci, refReads, altReads, iter, lhood, noisyErrorRate, cleanErrorRate, noisyLocusRate, simpleErrorRate, theta, \n";

    SingleSampleContextDataExportFormat exportedContextData;
    for (const auto& contextInfo : counts.getIndelCounts())
    {
        const auto& context(contextInfo.first);
        const auto& data(contextInfo.second);

        data.exportData(exportedContextData);

        if (exportedContextData.data.empty()) continue;

        log_os << "INFO: computing rates for context: " << context << "\n";
        reportExtendedContext(isLockTheta, context, exportedContextData, data, ros);
    }
}

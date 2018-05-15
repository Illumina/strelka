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

///
/// \author Chris Saunders
///

#include "snvModelVariantAndBinomialMixtureError.hh"

#include "blt_util/log.hh"
#include "blt_util/logSumUtil.hh"
#include "blt_util/prob_util.hh"
#include "blt_util/qscore.hh"

#define CODEMIN_USE_BOOST
#include "minimize_conj_direction.h"

#include "boost/math/special_functions/binomial.hpp"

#include <iomanip>
#include <iostream>
#include <memory>


//#define DEBUG_SNVMODEL3


using namespace BasecallCounts;


//static const double logCleanLocusBaseErrorRate(std::log(1e-7));
static
double
getLogCleanLocusBaseErrorRate(
    const double logNoisyRate,
    const double cleanErrorFactor)
{
    return logNoisyRate*cleanErrorFactor;
}



static
double
getObsLogLhood(
    const double logHomPrior,
    const double logHetPrior,
    const double logNoVariantPrior,
    const double* const logBaseErrorRate,
    const double* const logCleanBaseErrorRateByQual,
    const double logNoVariantRefRate,
    const double logNoVariantCleanRefRate,
    const double logNoisyLocusRate,
    const double logCleanLocusRate,
    const unsigned qualCount,
    const SingleSampleContextObservationPatternExportFormat& contextObservationPattern)
{
    static const double homAltRate(0.99);
    static const double hetAltRate(0.5);

    static const double logHomAltRate(std::log(homAltRate));
    static const double logHomRefRate(std::log1p(-homAltRate));
    static const double logHetRate(std::log(hetAltRate));

    const auto& s0(contextObservationPattern.strand0);
    const auto& s1(contextObservationPattern.strand1);

    const unsigned refQualTotal(s0.refAlleleCount+s1.refAlleleCount);

    // get lhood of homref GT:
    double noVariant_noise_s0(logNoVariantRefRate*s0.refAlleleCount);
    double noVariant_noise_s1(logNoVariantRefRate*s1.refAlleleCount);
    double noVariant_clean_s0(logNoVariantCleanRefRate*s0.refAlleleCount);
    double noVariant_clean_s1(logNoVariantCleanRefRate*s1.refAlleleCount);
    for (unsigned qualIndex(0); qualIndex<qualCount; ++qualIndex)
    {
        const double logCleanLocusBaseErrorRate(logCleanBaseErrorRateByQual[qualIndex]);

        noVariant_noise_s0 += logBaseErrorRate[qualIndex] * s0.altAlleleCount[qualIndex];
        noVariant_noise_s1 += logBaseErrorRate[qualIndex] * s1.altAlleleCount[qualIndex];
        noVariant_clean_s0 += logCleanLocusBaseErrorRate * s0.altAlleleCount[qualIndex];
        noVariant_clean_s1 += logCleanLocusBaseErrorRate * s1.altAlleleCount[qualIndex];
    }

    // first version treats the noisy locus as each strand, second version ignores strand and treats noisy locus by site:
#if 0
    // note we can typically ignore the multinomial coefficient because it is the same in all terms, but once we divide
    // by strand this term needs to be made explicit and accounted for:
    //
    // for practical reasons, see if we can get away with a binomial coefficient first, which ignores the differences between
    // different alt qualities:
    //

    // normally this usage of static wouldn't fly, but we're just trying this out:
    static bool isBcsComputed(false);
    static double logBcsnorm(0);

    if (! isBcsComputed)
    {
        unsigned altTotalStrand0(0);
        unsigned altTotalStrand1(0);
        for (unsigned qualIndex(0); qualIndex<qualCount; ++qualIndex)
        {
            altTotalStrand0 += s0.altCount[qualIndex];
            altTotalStrand1 += s1.altCount[qualIndex];
        }
        using namespace boost::math;
        const double bcs0(binomial_coefficient<double>((s0.refCount+altTotalStrand0),altTotalStrand0));
        const double bcs1(binomial_coefficient<double>((s1.refCount+altTotalStrand1),altTotalStrand1));
        const unsigned altQualTotal(altTotalStrand0+altTotalStrand1);
        const double bcsb(binomial_coefficient<double>((refQualTotal+altQualTotal), altQualTotal));
        logBcsnorm = std::log((bcs0*bcs1)/bcsb);
        isBcsComputed = true;
    }

    const double noVariant_s0(log_sum(logNoisyLocusRate + noVariant_noise_s0, logCleanLocusRate + noVariant_clean_s0));
    const double noVariant_s1(log_sum(logNoisyLocusRate + noVariant_noise_s1, logCleanLocusRate + noVariant_clean_s1));
    const double noVariant(noVariant_s0 + noVariant_s1 + logBcsnorm);
#else
    const double noVariant_noise(noVariant_noise_s0 + noVariant_noise_s1);
    const double noVariant_clean(noVariant_clean_s0 + noVariant_clean_s1);
    const double noVariant(getLogSum(logNoisyLocusRate + noVariant_noise, logCleanLocusRate + noVariant_clean));
#endif
    // get lhood of het GT:
    unsigned altQualTotal(0);
    for (unsigned qualIndex(0); qualIndex<qualCount; ++qualIndex)
    {
        altQualTotal += (s0.altAlleleCount[qualIndex] + s1.altAlleleCount[qualIndex]);
    }

    const double het(logHetRate*(refQualTotal+altQualTotal));

    // get lhood of hom GT:
    const double hom(logHomAltRate*altQualTotal + logHomRefRate*refQualTotal);

    return getLogSum(logHomPrior+hom, logHetPrior+het, logNoVariantPrior+noVariant);
}




static
double
contextLogLhood(
    const SingleSampleContextDataExportFormat& data,
    const double* const logBaseErrorRate,
    const double* const logCleanBaseErrorRateByQual,
    const double logNoisyLocusRate,
    const double logTheta)
{
    const unsigned qualCount(data.altAlleleBasecallErrorPhredProbLevels.size());

#ifdef DEBUG_SNVMODEL3
    log_os << "START loghood\n";
    log_os << "logBaseErrorRate:\n";
    for (unsigned qualIndex(0); qualIndex<qualCount; ++qualIndex)
    {
        log_os << logBaseErrorRate[qualIndex] << "\n";
    }
    log_os << "logCleanBaseErrorRateByQual:\n";
    for (unsigned qualIndex(0); qualIndex<qualCount; ++qualIndex)
    {
        log_os << logCleanBaseErrorRateByQual[qualIndex] << "\n";
    }
    log_os << "logNoisyLocusRate: " << logNoisyLocusRate << "\n";
    log_os << "logTheta: " << logTheta << "\n";
#endif

    static const double log2(std::log(2));
    const double logHomPrior(logTheta-log2);
    const double logHetPrior(logTheta);
    const double theta(std::exp(logTheta));
    const double logNoVariantPrior(std::log(1-(theta*3./2.)));


    uint64_t refTotal(0);
    double refErrorRateFactor(0);
    double cleanRefErrorRateFactor(0);
    for (unsigned qualIndex(0); qualIndex<qualCount; ++qualIndex)
    {
        const double logCleanLocusBaseErrorRate(logCleanBaseErrorRateByQual[qualIndex]);

        refErrorRateFactor += data.refCount[qualIndex] * std::exp(logBaseErrorRate[qualIndex]);
        cleanRefErrorRateFactor += data.refCount[qualIndex] * std::exp(logCleanLocusBaseErrorRate);
        refTotal += data.refCount[qualIndex];
    }
    const double noVariantRefRateComp(-(refErrorRateFactor / refTotal));
    const double logNoVariantRefRate(std::log1p(noVariantRefRateComp));

    const double noVariantCleanRefRateComp(-(cleanRefErrorRateFactor / refTotal));
    const double logNoVariantCleanRefRate(std::log1p(noVariantCleanRefRateComp));

    const double logCleanLocusRate(std::log1p(-std::exp(logNoisyLocusRate)));

    double logLhood(0.);
    for (const auto& value : data.observations)
    {
        const SingleSampleContextObservationPatternExportFormat& key(value.first);
        const unsigned repeatCount(value.second);

        const double mix(getObsLogLhood(
                             logHomPrior, logHetPrior, logNoVariantPrior, logBaseErrorRate, logCleanBaseErrorRateByQual, logNoVariantRefRate, logNoVariantCleanRefRate,
                             logNoisyLocusRate, logCleanLocusRate, qualCount, key));
        //const double cleanMix(getObsLogLhood(logHomPrior, logHetPrior, logAltHetPrior, logNoIndelPrior,
        //                                     logCleanLocusIndelRate, logCleanLocusIndelRate, logCleanLocusRefRate, obs));

        //const double mix(log_sum(logCleanLocusRate+cleanMix,logNoisyLocusRate+noisyMix));

        logLhood += (mix*repeatCount);
    }

#ifdef DEBUG_SNVMODEL3
    log_os << "EXIT loghood: " << logLhood << "\n";
#endif

    return logLhood;
}


// anon namespace for file-scoped structs:
namespace
{

struct error_minfunc : public codemin::minfunc_interface<double>
{
    explicit
    error_minfunc(
        const SingleSampleContextDataExportFormat& data,
        const bool isFreeCleanLocusError,
        const bool isLockTheta)
        : _data(data),
          _isFreeCleanLocusError(isFreeCleanLocusError),
          _isLockTheta(isLockTheta),
          _qualParamSize(_data.altAlleleBasecallErrorPhredProbLevels.size()),
          _args(new double[getArgCount()])
    {}

    virtual unsigned dim() const
    {
        return getParamCount();
    }

    virtual double val(const double* in)
    {
        paramsToArgs(in,_args.get());
        return -contextLogLhood(_data,
                                _args.get(),
                                _args.get()+_qualParamSize,
                                _args[getNoisyLocusArgIndex()],
                                _args[getThetaArgIndex()]);
    }

    /// normalize the minimization values back to usable function arguments
    ///
    /// most values are not valid on [-inf,inf] -- the minimizer doesn't
    /// know this. here is where we fill in the gap:
    ///
    void
    paramsToArgs(
        const double* in,
        double* out) const
    {
        auto rateSmoother = [](double a) -> double
        {
            static const double triggerVal(1e-2);
            static const double limitVal(0.5);
            static const double logTriggerVal(std::log(triggerVal));
            static const double logLimitVal(std::log(limitVal));
            if (a>logTriggerVal)
            {
                a = std::log1p(a-logTriggerVal) + logTriggerVal;
            }
            return (a>logLimitVal ? logLimitVal-std::abs(a-logLimitVal) : a);
        };

        // A lot of conditioning is required to keep the model from winding
        // theta around zero and getting confused, here we start applying a
        // second log to the delta above triggerTheta, and finally put a hard stop
        // at maxTheta -- hard stops are obviously bad b/c the model can get lost
        // on the flat plane even if the ML value is well below this limit, but
        // in practice this is such a ridiculously high value for theta, that
        // I don't see the model getting trapped.
        auto thetaSmoother = [](double a) -> double
        {
            static const double triggerVal(1e-2);
            static const double limitVal(0.3);
            static const double logTriggerVal(std::log(triggerVal));
            static const double logLimitVal(std::log(limitVal));
            if (a>logTriggerVal)
            {
                a = std::log1p(a-logTriggerVal) + logTriggerVal;
            }
            return (a>logLimitVal ? logLimitVal-std::abs(a-logLimitVal) : a);
        };

        if (_isFreeCleanLocusError)
        {
            const unsigned qualParamCount(getQualParamCount());
            for (unsigned qualIndex(0); qualIndex<qualParamCount; ++qualIndex)
            {
                out[qualIndex] =  rateSmoother(in[qualIndex]);
            }
        }
        else
        {
            for (unsigned qualIndex(0); qualIndex<_qualParamSize; ++qualIndex)
            {
                out[qualIndex] =  rateSmoother(in[qualIndex]);
            }
            for (unsigned qualIndex(0); qualIndex<_qualParamSize; ++qualIndex)
            {
                out[qualIndex+_qualParamSize] = getLogCleanLocusBaseErrorRate(out[qualIndex], in[_qualParamSize]);
            }
        }
        out[getNoisyLocusArgIndex()] = rateSmoother(in[getNoisyLocusParamIndex()]);

        if (_isLockTheta)
        {
            out[getThetaArgIndex()] = defaultLogTheta;
        }
        else
        {
            out[getThetaArgIndex()] = thetaSmoother(in[getThetaParamIndex()]);
        }
    }

    unsigned
    getQualParamCount() const
    {
        if (_isFreeCleanLocusError)
        {
            return _qualParamSize*2;
        }
        else
        {
            return _qualParamSize+1;
        }
    }

    unsigned
    getQualArgCount() const
    {
        return _qualParamSize*2;
    }

    unsigned
    getNoisyLocusParamIndex() const
    {
        return getQualParamCount();
    }

    unsigned
    getThetaParamIndex() const
    {
        return getQualParamCount()+1;
    }

    unsigned
    getNoisyLocusArgIndex() const
    {
        return getQualArgCount();
    }

    unsigned
    getThetaArgIndex() const
    {
        return getQualArgCount()+1;
    }

    unsigned
    getParamCount() const
    {
        return (getQualParamCount() + 1 + (_isLockTheta ? 0 : 1));
    }

    unsigned
    getArgCount() const
    {
        return (getQualArgCount() + 2);
    }

    static const double defaultLogTheta;

private:
    const SingleSampleContextDataExportFormat& _data;
    bool _isFreeCleanLocusError;
    bool _isLockTheta;
    unsigned _qualParamSize;
    std::unique_ptr<double[]> _args;
};

const double error_minfunc::defaultLogTheta = std::log(1e-3);


struct SignalGroupTotal
{
    double locus = 0;
    double skipped = 0;
    double ref = 0;
    double alt = 0;
};

} // END anon namespace for file-scoped structs



static
void
reportQualErrorRateSet(
    const Context& context,
    const ContextData& contextData,
    const uint16_t qual,
    const SignalGroupTotal& sigTotal,
    unsigned iter,
    const double loghood,
    const double noisyErrorRate,
    const double cleanErrorRate,
    const double noisyLocusRate,
    const double theta,
    const bool isFreeCleanLocusError,
    const double cleanErrorFactor,
    std::ostream& os)
{
    static const std::string sep(", ");

    const double expectErrorRate(qphred_to_error_prob(qual));

    os << std::setprecision(10);
    os << context
       << sep << "Q" << qual
       << sep << contextData.excludedRegionSkipped
       << sep << (sigTotal.skipped+sigTotal.locus)
       << sep << sigTotal.locus
       << sep << sigTotal.ref
       << sep << sigTotal.alt
       << sep << iter
       << sep << loghood
       << sep << noisyErrorRate
       << sep << cleanErrorRate
       << sep;

    if (isFreeCleanLocusError)
    {
        os << "NA";
    }
    else
    {
        os << cleanErrorFactor;
    }

    os << sep << noisyLocusRate
       << sep << (noisyErrorRate*noisyLocusRate + cleanErrorRate*(1-noisyLocusRate))
       << sep << expectErrorRate
       << sep << theta
       << "\n";
}



static
void
reportExtendedContext(
    const bool isFreeCleanLocusError,
    const bool isLockTheta,
    const Context& context,
    const ContextData& contextData,
    std::ostream& os)
{
    SingleSampleContextDataExportFormat exportData;
    contextData.counts.exportData(exportData);
    error_minfunc errFunc(exportData,isFreeCleanLocusError,isLockTheta);

    const unsigned qualCount(exportData.altAlleleBasecallErrorPhredProbLevels.size());
    const unsigned paramCount(errFunc.getParamCount());
    const unsigned argCount(errFunc.getArgCount());

    // initialize conjugate direction minimizer settings and minimize lhood...
    //
    std::unique_ptr<double[]> minParams(new double[paramCount]);
    std::unique_ptr<double[]> minArgs(new double[argCount]);

    unsigned iter;
    double x_all_loghood;
    {
        static const double line_tol(1e-10);
        static const double end_tol(1e-10);
        static const unsigned max_iter(40);

        // initialize parameter search
        for (unsigned qualIndex(0); qualIndex<qualCount; ++qualIndex)
        {
            minParams[qualIndex] = std::log(1e-2);
        }
        if (isFreeCleanLocusError)
        {
            for (unsigned qualIndex(qualCount); qualIndex<qualCount*2; ++qualIndex)
            {
                minParams[qualIndex] = std::log(1e-4);
            }
        }
        else
        {
            minParams[qualCount] = 2;
        }
        minParams[errFunc.getNoisyLocusParamIndex()] = std::log(1e-2);
        if (! isLockTheta)
        {
            minParams[errFunc.getThetaParamIndex()] = error_minfunc::defaultLogTheta;
        }

        const unsigned paramCount2(paramCount*paramCount);
        std::unique_ptr<double[]> conjDir(new double[paramCount2]);

        std::fill(conjDir.get(),conjDir.get()+paramCount2,0.);
        const unsigned dim(errFunc.dim());
        for (unsigned i(0); i<dim; ++i)
        {
            conjDir[i*(dim+1)] = 0.001;
        }

        double start_tol(end_tol);
        double final_dlh;

        codemin::minimize_conj_direction(minParams.get(),conjDir.get(),errFunc,start_tol,end_tol,line_tol,
                                         x_all_loghood,iter,final_dlh,max_iter);

        errFunc.paramsToArgs(minParams.get(),minArgs.get());
    }

    // report:
    {
        const double noisyLocusRate(std::exp(minArgs[errFunc.getNoisyLocusArgIndex()]));
        const double theta(std::exp(minArgs[errFunc.getThetaArgIndex()]));

        SignalGroupTotal sigTotal;
        {
            const uint64_t emptySkippedLociCount(contextData.emptySkipped);
            const uint64_t depthSkippedLociCount(contextData.depthSkipped);
            const uint64_t noisyContextSkippedLociCount(contextData.noiseSkipped);

            sigTotal.skipped = (emptySkippedLociCount+depthSkippedLociCount+noisyContextSkippedLociCount);

            for (const auto& value : exportData.observations)
            {
                const unsigned repeatCount(value.second);
                sigTotal.locus += repeatCount;
            }
        }

        for (unsigned qualIndex(0); qualIndex<qualCount; ++qualIndex)
        {
            sigTotal.ref = exportData.refCount[qualIndex];
            sigTotal.alt = 0;
            for (const auto& value : exportData.observations)
            {
                const auto& key(value.first);
                const unsigned repeatCount(value.second);
                const unsigned altCount(key.strand0.altAlleleCount[qualIndex]+key.strand1.altAlleleCount[qualIndex]);
                sigTotal.alt += (altCount*repeatCount);
            }

            const double noisyErrorRate(std::exp(minArgs[qualIndex]));
//            const double logCleanLocusBaseErrorRate(getLogCleanLocusBaseErrorRate(normParams[qualIndex]));

            double cleanErrorFactor(0);
            double cleanErrorRate(0);
            if (isFreeCleanLocusError)
            {
                cleanErrorRate = std::exp(minArgs[qualIndex+qualCount]);
            }
            else
            {
                cleanErrorFactor = minParams[qualCount];
                cleanErrorRate = std::exp(getLogCleanLocusBaseErrorRate(minArgs[qualIndex],cleanErrorFactor));
            }
            reportQualErrorRateSet(context, contextData, exportData.altAlleleBasecallErrorPhredProbLevels[qualIndex], sigTotal, iter, -x_all_loghood,
                                   noisyErrorRate, cleanErrorRate, noisyLocusRate, theta, isFreeCleanLocusError, cleanErrorFactor, os);
        }
    }
}



void
snvModelVariantAndBinomialMixtureError(
    const SequenceAlleleCounts& counts)
{
    const bool isFreeCleanLocusError(false);
    const bool isLockTheta(true);

    std::ostream& ros(std::cout);

    ros << "context, qual, excludedLoci, nonExcludedLoci, usedLoci, refReads, altReads, iter, lhood, noisyErrorRate, cleanErrorRate, cleanErrorFactor, noisyLocusRate, simpleErrorRate, expectErrorRate, theta\n";

    for (const auto& contextInfo : counts.getBaseCounts())
    {
        const auto& context(contextInfo.first);
        const auto& data(contextInfo.second);

        reportExtendedContext(isFreeCleanLocusError,isLockTheta, context, data, ros);
    }
}

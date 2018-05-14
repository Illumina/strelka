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

#include "snvModelVariantAndIndyError.hh"

#include "blt_util/logSumUtil.hh"
#include "blt_util/prob_util.hh"
#include "blt_util/qscore.hh"

#define CODEMIN_USE_BOOST
#include "minimize_conj_direction.h"

#include <iomanip>
#include <iostream>
#include <memory>


using namespace BasecallCounts;


static
double
contextLogLhood(
    const SingleSampleContextDataExportFormat& data,
    const double* logBaseErrorRate,
    const double logTheta)
{
    static const double homAltRate(0.99);
    static const double hetAltRate(0.5);

    static const double logHomAltRate(std::log(homAltRate));
    static const double logHomRefRate(std::log(1.-homAltRate));
    static const double logHetRate(std::log(hetAltRate));

    static const double log2(std::log(2));
    const double logHomPrior(logTheta-log2);
    const double logHetPrior(logTheta);
    const double theta(std::exp(logTheta));
    const double logNoVariantPrior(std::log(1-(theta*3./2.)));

    const unsigned qualCount(data.altAlleleBasecallErrorPhredProbLevels.size());

    uint64_t refTotal(0);
    double refErrorRateFactor(0);
    for (unsigned qualIndex(0); qualIndex<qualCount; ++qualIndex)
    {
        refErrorRateFactor += data.refCount[qualIndex] * std::exp(logBaseErrorRate[qualIndex]);
        refTotal += data.refCount[qualIndex];
    }
    const double noVariantRefRate(1. - (refErrorRateFactor / refTotal));
    const double logNoVariantRefRate(std::log(noVariantRefRate));

    double logLhood(0.);
    for (const auto& value : data.observations)
    {
        const auto& key(value.first);
        const unsigned repeatCount(value.second);

        const auto& s0(key.strand0);
        const auto& s1(key.strand1);

        const unsigned refQualTotal(s0.refAlleleCount+s1.refAlleleCount);

        // get lhood of homref GT:
        double noVariant(logNoVariantRefRate*refQualTotal);
        for (unsigned qualIndex(0); qualIndex<qualCount; ++qualIndex)
        {
            noVariant += logBaseErrorRate[qualIndex] * (s0.altAlleleCount[qualIndex] + s1.altAlleleCount[qualIndex]);
        }

        // get lhood of het GT:
        unsigned altQualTotal(0);
        for (unsigned qualIndex(0); qualIndex<qualCount; ++qualIndex)
        {
            altQualTotal += (s0.altAlleleCount[qualIndex] + s1.altAlleleCount[qualIndex]);
        }

        const double het(logHetRate*(refQualTotal+altQualTotal));

        // get lhood of hom GT:
        const double hom(logHomAltRate*altQualTotal + logHomRefRate*refQualTotal);

        const double mix(getLogSum(logHomPrior+hom, logHetPrior+het, logNoVariantPrior+noVariant));

        logLhood += (mix*repeatCount);
    }

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
        const bool isLockTheta = false)
        : _data(data),
          _isLockTheta(isLockTheta),
          _qualParamSize(_data.altAlleleBasecallErrorPhredProbLevels.size()),
          _params(new double[_qualParamSize+1])
    {}

    virtual unsigned dim() const
    {
        return _qualParamSize + (_isLockTheta ? 0 : 1);
    }

    virtual double val(const double* in)
    {
        argToParameters(in,_params.get());
        return -contextLogLhood(_data,
                                _params.get(),
                                (_isLockTheta ? defaultLogTheta : _params[_qualParamSize]));
    }

    /// normalize the minimization values back to usable parameters
    ///
    /// most values are not valid on [-inf,inf] -- the minimizer doesn't
    /// know this. here is where we fill in the gap:
    ///
    void
    argToParameters(
        const double* in,
        double* out)
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

        for (unsigned qualIndex(0); qualIndex<_qualParamSize; ++qualIndex)
        {
            out[qualIndex] =  rateSmoother(in[qualIndex]);
        }
        out[_qualParamSize] = thetaSmoother(in[_qualParamSize]);
    }

    static const double defaultLogTheta;

private:
    const SingleSampleContextDataExportFormat& _data;
    bool _isLockTheta;
    unsigned _qualParamSize;
    std::unique_ptr<double[]> _params;
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
    const double errorRate,
    const double theta,
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
       << sep << errorRate
       << sep << expectErrorRate
       << sep << theta
       << "\n";
}



static
void
reportExtendedContext(
    const bool isLockTheta,
    const Context& context,
    const ContextData& contextData,
    std::ostream& os)
{
    SingleSampleContextDataExportFormat exportData;
    contextData.counts.exportData(exportData);

    const unsigned qualCount(exportData.altAlleleBasecallErrorPhredProbLevels.size());
    const unsigned paramCount(qualCount+1);

    // initialize conjugate direction minimizer settings and minimize lhood...
    //
    std::unique_ptr<double[]> minParams(new double[paramCount]);
    std::unique_ptr<double[]> normParams(new double[paramCount]);

    unsigned iter;
    double x_all_loghood;
    {
        static const double line_tol(1e-10);
        static const double end_tol(1e-10);
        static const unsigned max_iter(20);

        // initialize parameter search
        for (unsigned qualIndex(0); qualIndex<qualCount; ++qualIndex)
        {
            minParams[qualIndex] = std::log(1e-4);
        }
        minParams[qualCount] = error_minfunc::defaultLogTheta;

        const unsigned paramCount2(paramCount*paramCount);
        std::unique_ptr<double[]> conjDir(new double[paramCount2]);

        std::fill(conjDir.get(),conjDir.get()+paramCount2,0.);
        const unsigned dim(qualCount + (isLockTheta ? 0 : 1));
        for (unsigned i(0); i<dim; ++i)
        {
            conjDir[i*(dim+1)] = 0.001;
        }

        double start_tol(end_tol);
        double final_dlh;
        error_minfunc errFunc(exportData,isLockTheta);

        codemin::minimize_conj_direction(minParams.get(),conjDir.get(),errFunc,start_tol,end_tol,line_tol,
                                         x_all_loghood,iter,final_dlh,max_iter);

        errFunc.argToParameters(minParams.get(),normParams.get());
    }

    // report:
    {
        const double theta(std::exp(normParams[qualCount]));

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

            const double qualErrorRate(std::exp(normParams[qualIndex]));
            reportQualErrorRateSet(context, contextData, exportData.altAlleleBasecallErrorPhredProbLevels[qualIndex], sigTotal, iter, -x_all_loghood, qualErrorRate, theta, os);
        }
    }
}



void
snvModelVariantAndIndyError(
    const SequenceAlleleCounts& counts)
{
    const bool isLockTheta(false);

    std::ostream& ros(std::cout);

    ros << "context, qual, excludedLoci, nonExcludedLoci, usedLoci, refReads, altReads, iter, lhood, errorRate, expectErrorRate, theta\n";

    for (const auto& contextInfo : counts.getBaseCounts())
    {
        const auto& context(contextInfo.first);
        const auto& data(contextInfo.second);

        reportExtendedContext(isLockTheta, context, data, ros);
    }
}

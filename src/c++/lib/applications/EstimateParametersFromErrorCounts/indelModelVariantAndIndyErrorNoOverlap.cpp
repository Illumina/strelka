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

#include "indelModelVariantAndIndyErrorNoOverlap.hh"

#include "blt_util/logSumUtil.hh"
#include "blt_util/prob_util.hh"

#define CODEMIN_USE_BOOST
#include "minimize_conj_direction.h"

#include <iomanip>
#include <iostream>

namespace
{

namespace MIN_PARAMS
{
enum index_t
{
    LN_INDEL_ERROR_RATE,
    LN_THETA,
    SIZE
};
}


using namespace IndelCounts;



static
double
contextLogLhood(
    const SingleSampleContextDataExportFormat& exportedContextData,
    const double logIndelErrorRate,
    const bool isInsert,
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
    const double logNoIndelPrior(std::log(1-(theta*3./2.)));

    const double logNoIndelRefRate(std::log1p(-std::exp(logIndelErrorRate)));

    double logLhood(0.);
    for (const auto& contextObservationInfo : exportedContextData.data)
    {
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

        const double hom(logHomAltRate* +totalIndelObservations +
                         logHomRefRate*contextObservationInfo.refObservations);

        const double mix(getLogSum(logHomPrior+hom, logHetPrior+het, logNoIndelPrior+noindel));

        logLhood += (mix*contextObservationInfo.contextInstanceCount);
    }

    return logLhood;
}


struct error_minfunc : public codemin::minfunc_interface<double>
{
    explicit
    error_minfunc(
        const SingleSampleContextDataExportFormat& exportedContextData,
        const bool isInsert,
        const bool isLockTheta = false)
        : _exportedContextData(exportedContextData),
          _isInsert(isInsert),
          _isLockTheta(isLockTheta)
    {}

    virtual unsigned dim() const
    {
        return (_isLockTheta ? (MIN_PARAMS::SIZE-1) : MIN_PARAMS::SIZE);
    }

    virtual double val(const double* in)
    {
        // std::cerr << "Submitting: " << in[0] << " " << in[1] << " " << in[2] << "\n";
        argToParameters(in,_params);
        //  std::cerr << "Trying: ins/del/theta: " << std::exp(_params[0]) << " " << std::exp(_params[1]) << " " << std::exp(_params[2]) << "\n";
        return -contextLogLhood(_exportedContextData,
                                _params[MIN_PARAMS::LN_INDEL_ERROR_RATE],
                                _isInsert,
                                (_isLockTheta ? defaultLogTheta : _params[MIN_PARAMS::LN_THETA]));
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
            static const double triggerVal(1e-3);
            static const double limitVal(0.3);
            static const double logTriggerVal(std::log(triggerVal));
            static const double logLimitVal(std::log(limitVal));
            if (a>logTriggerVal)
            {
                a = std::log1p(a-logTriggerVal) + logTriggerVal;
            }
            return (a>logLimitVal ? logLimitVal-std::abs(a-logLimitVal) : a);
        };

        out[MIN_PARAMS::LN_INDEL_ERROR_RATE] =  rateSmoother(in[MIN_PARAMS::LN_INDEL_ERROR_RATE]);
        out[MIN_PARAMS::LN_THETA] = thetaSmoother(in[MIN_PARAMS::LN_THETA]);
    }

    static const double defaultLogTheta;

private:
    const SingleSampleContextDataExportFormat& _exportedContextData;
    bool _isInsert;
    bool _isLockTheta;
    double _params[MIN_PARAMS::SIZE];
};

const double error_minfunc::defaultLogTheta = std::log(1e-4);


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
    const double indelErrorRate,
    const double theta,
    std::ostream& os)
{
    static const std::string sep(", ");

    os << std::setprecision(10);
    os << context << "_" << extendedContextTag << sep
       << contextData.excludedRegionSkipped << sep
       << (sigTotal.locus + contextData.depthSkipped) << sep
       << sigTotal.locus << sep
       << sigTotal.ref << sep
       << sigTotal.alt << sep
       << iter << sep
       << loghood << sep
       << indelErrorRate << sep
       << theta << "\n";
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

        double minParams[MIN_PARAMS::SIZE];

        unsigned iter;
        double x_all_loghood;
        {
            static const double line_tol(1e-10);
            static const double end_tol(1e-10);
            static const unsigned max_iter(20);

            // initialize parameter search
            minParams[MIN_PARAMS::LN_INDEL_ERROR_RATE] = std::log(1e-3);
            minParams[MIN_PARAMS::LN_THETA] = error_minfunc::defaultLogTheta;

            static const unsigned SIZE2(MIN_PARAMS::SIZE*MIN_PARAMS::SIZE);
            double conjDir[SIZE2];

            std::fill(conjDir,conjDir+SIZE2,0.);
            const unsigned dim(isLockTheta ? MIN_PARAMS::SIZE-1 : MIN_PARAMS::SIZE);
            for (unsigned dimIndex(0); dimIndex<dim; ++dimIndex)
            {
                conjDir[dimIndex*(dim+1)] = 0.001;
            }

            double start_tol(end_tol);
            double final_dlh;
            error_minfunc errFunc(exportedContextData, isInsert, isLockTheta);

            codemin::minimize_conj_direction(minParams,conjDir,errFunc,start_tol,end_tol,line_tol,
                                             x_all_loghood,iter,final_dlh,max_iter);
        }

        // report:
        {
            double normalizedParams[MIN_PARAMS::SIZE];
            error_minfunc::argToParameters(minParams,normalizedParams);

            const double theta(std::exp(normalizedParams[MIN_PARAMS::LN_THETA]));

            const double indelErrorRate(std::exp(normalizedParams[MIN_PARAMS::LN_INDEL_ERROR_RATE]));
            const std::string tag(isInsert ? "I" : "D");
            reportIndelErrorRateSet(context, tag.c_str(), sigInsertTotal, contextData, iter, -x_all_loghood, indelErrorRate, theta, os);
        }
    }
}

}



void
indelModelVariantAndIndyErrorNoOverlap(
    const SequenceAlleleCounts& counts)
{
    const bool isLockTheta(false);

    std::ostream& ros(std::cout);

    ros << "context, excludedLoci, nonExcludedLoci, usedLoci, refReads, altReads, iter, lhood, rate, theta\n";

    SingleSampleContextDataExportFormat exportedContextData;
    for (const auto& contextInfo : counts.getIndelCounts())
    {
        const auto& context(contextInfo.first);
        const auto& data(contextInfo.second);

        data.exportData(exportedContextData);

        if (exportedContextData.data.empty()) continue;

        std::cerr << "INFO: computing rates for context: " << context << "\n";
        reportExtendedContext(isLockTheta, context, exportedContextData, data, ros);
    }
}

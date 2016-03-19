// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2016 Illumina, Inc.
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

#include "modelVariantAndIndyError.hh"

#include "blt_util/math_util.hh"
#include "blt_util/prob_util.hh"

#define CODEMIN_USE_BOOST
#include "minimize_conj_direction.h"

#include <iomanip>
#include <iostream>

namespace MIN_PARAMS
{
    enum index_t {
        LN_ALT_ERROR_RATE,
        LN_THETA,
        SIZE
    };
}



static
double
contextLogLhood(
    const std::vector<ExportedObservations>& observations,
    const double logNoIndelAltRate,
    const double logTheta)
{
    static const double homAltRate(0.99);
    static const double hetAltRate(0.5);

    static const double logHomAltRate(std::log(homAltRate));
    static const double logHomRefRate(std::log(1.-homAltRate));
    static const double logHetRate(std::log(hetAltRate));

    static const double ln2(std::log(2));
    const double logHomPrior(logTheta-ln2);
    const double logHetPrior(logTheta);
    const double logNoIndelPrior(std::log(1-std::exp(logTheta)*3./2.));

    const double logNoIndelRefRate(std::log(1-std::exp(logNoIndelAltRate)));

    double logLhood(0.);
    for (const auto& obs : observations)
    {
        const unsigned totalObservations(obs.altObservations+obs.refObservations);
        const double het(logHetRate*totalObservations);
        const double hom(logHomAltRate*obs.altObservations + logHomRefRate*obs.refObservations);
        const double noindel(logNoIndelAltRate*obs.altObservations + logNoIndelRefRate*obs.refObservations);

        /// TODO: generalize log_sum to N values...
        const double mix(log_sum(log_sum(logHomPrior+hom,logHetPrior+het),logNoIndelPrior+noindel));

        logLhood += (mix*obs.repeatCount);
    }

    return logLhood;
}


struct error_minfunc : public codemin::minfunc_interface<double>
{
    explicit
    error_minfunc(
        const std::vector<ExportedObservations>& observations,
        const bool isLockTheta = false)
    : _obs(observations), _isLockTheta(isLockTheta)
    {}

    virtual unsigned dim() const
    {
        return (_isLockTheta ? 1 : MIN_PARAMS::SIZE);
    }

    virtual double val(const double* in)
    {
        argToParameters(in,_params);
        return -contextLogLhood(_obs,
                _params[MIN_PARAMS::LN_ALT_ERROR_RATE],
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
        // keep any log(prob) value negative in case the minimizer bends it
        // around the corner
        auto convert = [](const double a) -> double { return -std::abs(a); };

        // A lot of conditioning is required to keep the model from winding
        // theta around zero and getting confused, here we start applying a
        // second log to the delta above 1e-2, and finally put a hard stop
        // at 0.4 -- hard stops are obviously bad b/c the model can get lost
        // on the flat plane even if the ML value is well below this limit, but
        // in practice this is such a ridiculously high value for theta, that
        // I don't see the model getting trapped.
        auto thetaSmoother = [](double a) -> double {
            static const double limitVal(std::log(0.4));
            static const double triggerVal(std::log(1e-2));
            if (a>triggerVal)
            {
                a = std::log(1+(a-triggerVal)) + triggerVal;
            }
            return (a>limitVal ? limitVal : a);
        };

        out[MIN_PARAMS::LN_ALT_ERROR_RATE] =  convert(in[MIN_PARAMS::LN_ALT_ERROR_RATE]);
        out[MIN_PARAMS::LN_THETA] = thetaSmoother(in[MIN_PARAMS::LN_THETA]);
    }

    static constexpr double defaultLogTheta = std::log(1e-4);

private:
    const std::vector<ExportedObservations>& _obs;
    bool _isLockTheta;
    double _params[MIN_PARAMS::SIZE];
};


void
modelVariantAndIndyError(
    const SequenceErrorCounts& counts)
{
    const bool isLockTheta(false);

    std::ostream& ros(std::cout);

    ros << "context, loci, reads, iter, lhood, rate, theta\n";

    double minParams[MIN_PARAMS::SIZE];
    static const unsigned SIZE2(MIN_PARAMS::SIZE*MIN_PARAMS::SIZE);
    double conjDir[SIZE2];

    std::vector<ExportedObservations> observations;
    for (const auto& contextInfo : counts)
    {
        const auto& context(contextInfo.first);
        const auto& data(contextInfo.second);

        data.exportObservations(observations);

        // get some useful counts unrelated to minimization or model:
        double refTotal(0.);
        double altTotal(0.);
        double locusTotal(0.);

        for (const auto& obs : observations)
        {
            refTotal += (obs.refObservations*obs.repeatCount);
            altTotal += (obs.altObservations*obs.repeatCount);
            locusTotal += obs.repeatCount;
        }

        /// initialize conj dir minimizer settings and min...
        ///
        unsigned iter;
        double x_all_loghood;
        {
            static const double line_tol(1e-10);
            static const double end_tol(1e-10);
            static const unsigned max_iter(200);

            // initialize parameter search
            minParams[MIN_PARAMS::LN_ALT_ERROR_RATE] = std::log(1e-4);
            minParams[MIN_PARAMS::LN_THETA] = error_minfunc::defaultLogTheta;

            std::fill(conjDir,conjDir+SIZE2,0.);
            for (unsigned i(0); i<MIN_PARAMS::SIZE; ++i)
            {
                conjDir[i*(MIN_PARAMS::SIZE+1)] = 0.01;
            }

            double start_tol(end_tol);
            double final_dlh;
            error_minfunc errFunc(observations,isLockTheta);

            codemin::minimize_conj_direction(minParams,conjDir,errFunc,start_tol,end_tol,line_tol,
                                             x_all_loghood,iter,final_dlh,max_iter);
        }

        double normalizedParams[MIN_PARAMS::SIZE];
        error_minfunc::argToParameters(minParams,normalizedParams);

        const double NoIndelAltRate(std::exp(normalizedParams[MIN_PARAMS::LN_ALT_ERROR_RATE]));
        const double theta(std::exp(normalizedParams[MIN_PARAMS::LN_THETA]));

        {
            static const std::string sep(", ");

            ros << std::setprecision(10);
            ros << context << sep
                << locusTotal << sep
                << (altTotal+refTotal) << sep
                << iter << sep
                << -x_all_loghood << sep
                << NoIndelAltRate << sep
                << theta << "\n";
                ;
        }
    }
}

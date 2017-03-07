// -*- mode: c++; indent-tabs-mode: nil; -*-
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

#include "starling_continuous_variant_caller.hh"
#include "blt_util/qscore.hh"

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/distributions/binomial.hpp>



static double AssignPValue(unsigned observedCallCount, unsigned coverage, unsigned estimatedBaseCallQuality)
{
    if (observedCallCount == 0)
        return 1.0;

    const double errorRate = qphred_to_error_prob(estimatedBaseCallQuality);
    return (boost::math::gamma_p(observedCallCount, coverage * errorRate));
}



int
starling_continuous_variant_caller::
poisson_qscore(unsigned callCount, unsigned coverage, unsigned estimatedBaseCallQuality, int maxQScore)
{
    double pValue = AssignPValue(callCount, coverage, estimatedBaseCallQuality);
    if (pValue <= 0) return maxQScore;
    return std::min(maxQScore, error_prob_to_qphred(pValue));
}



static
double
binomialLogDensity(
    unsigned trials,
    unsigned successes,
    double successProb)
{
    using namespace boost::math;

    assert((successProb >= 0.) and (successProb <= 1.));
    assert(successes <= trials);

    if (trials==0) return 0;
    return std::log(pdf(binomial(trials, successProb), successes));
}



double
starling_continuous_variant_caller::
strandBias(
    unsigned fwdAlt,
    unsigned revAlt,
    unsigned fwdOther,
    unsigned revOther)
{
    const unsigned fwdTotal(fwdAlt+fwdOther);
    const unsigned revTotal(revAlt+revOther);
    const unsigned total(fwdTotal+revTotal);
    if (total==0) return 0;

    const double fwdAltFreq(safeFrac(fwdAlt,fwdTotal));
    const double revAltFreq(safeFrac(revAlt,revTotal));
    const double altFreq(safeFrac(fwdAlt + revAlt, total));

    static const double errorRate(0.005);

    const double fwdLnp(binomialLogDensity( fwdTotal, fwdAlt, fwdAltFreq) + binomialLogDensity( revTotal, revAlt, errorRate));
    const double revLnp(binomialLogDensity( fwdTotal, fwdAlt, errorRate) + binomialLogDensity( revTotal, revAlt, revAltFreq));
    const double lnp(binomialLogDensity( fwdTotal, fwdAlt, altFreq) + binomialLogDensity( revTotal, revAlt, altFreq));
    return std::max(-100., std::max(fwdLnp, revLnp) - lnp);
}

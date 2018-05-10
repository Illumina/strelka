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

#pragma once

#include "starling_base_shared.hh"
#include "blt_util/logSumUtil.hh"


/// Return the approximate log likelihood of the read being produced from a random location in the genome
inline
double
getIncorrectMappingLogLikelihood(
    const starling_base_deriv_options& dopt,
    const bool isTier2,
    const uint16_t nonAmbiguousBasesInRead)
{
    const double thisRandomBaseMatchLogLikelihood(isTier2 ?
                                                  dopt.tier2RandomBaseMatchLogProb :
                                                  dopt.randomBaseMatchLogProb );
    return thisRandomBaseMatchLogLikelihood*nonAmbiguousBasesInRead;
}


/// Given the log-likelihood of the read conditioned on correct mapping as input, sum over the
/// correct and incorrect mapping states to approximately integrate out the mapping status condition.
///
inline
double
integrateOutMappingStatus(
    const starling_base_deriv_options& dopt,
    const uint16_t nonAmbiguousBasesInRead,
    const double correctMappingLogLikelihood,
    const bool isTier2)
{
    // the second term formally has an incorrect mapping prior (prior of incorrectly mapping a read at random),
    // but this is effectively 1 so it is approximated out
    return getLogSum((correctMappingLogLikelihood + dopt.correctMappingLogPrior),
                     getIncorrectMappingLogLikelihood(dopt, isTier2, nonAmbiguousBasesInRead) /* + incorrectMappingLogPrior ~= 0 */);
}

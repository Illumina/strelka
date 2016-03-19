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

#include "model1.hh"

#include "blt_util/math_util.hh"

#include "boost/math/distributions/binomial.hpp"

#include <iomanip>
#include <iostream>



void
model1(
    const SequenceErrorCounts& counts)
{
    static const unsigned minDepth(25);
    static const double maxAltFrac(0.05);

    std::ostream& ros(std::cout);

    ros << "context, loci, reads, rate, rate_95%_upper_bound\n";

    std::vector<ExportedObservations> observations;
    for (const auto& contextInfo : counts)
    {
        const auto& context(contextInfo.first);
        const auto& data(contextInfo.second);

        data.exportObservations(observations);

        double refTotal(0.);
        double altTotal(0.);
        double locusTotal(0.);

        for (const auto& obs : observations)
        {
            const unsigned total(obs.refObservations+obs.altObservations);
            const double altFrac(static_cast<double>(obs.altObservations)/total);
            if ( (total<minDepth) || (altFrac > maxAltFrac)) continue;

            refTotal += (obs.refObservations*obs.repeatCount);
            altTotal += (obs.altObservations*obs.repeatCount);
            locusTotal += obs.repeatCount;
        }

        {
            static const std::string sep(", ");
            const double total(refTotal+altTotal);

            static const double alpha(0.05);

            const double upper(boost::math::binomial_distribution<double>::find_upper_bound_on_p(total, altTotal, alpha));

            ros << std::setprecision(10);
            ros << context << sep
                << locusTotal << sep
                << total << sep
                << safeFrac(altTotal,total) << sep
                << upper << "\n";
                ;
        }
    }
}

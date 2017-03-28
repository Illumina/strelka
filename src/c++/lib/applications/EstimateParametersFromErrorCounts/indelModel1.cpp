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

#include "blt_util/math_util.hh"

#include "boost/math/distributions/binomial.hpp"

#include <iomanip>
#include <iostream>
#include <numeric>
#include "indelModel1.hh"

namespace
{


struct SignalGroupTotal
{
    double ref = 0;
    double alt = 0;
    double noiseLocus = 0;
    double locus = 0;
};



static
void
reportExtendedContext(
    const double maxAltFrac,
    const unsigned minDepth,
    const IndelErrorContext& context,
    const std::vector<ExportedIndelObservations>& observations,
    const IndelErrorData& data,
    const unsigned altBeginIndex,
    const unsigned altEndIndex,
    const char* extendedContextTag,
    std::ostream& os)
{
    SignalGroupTotal sigTotal;

    for (const ExportedIndelObservations& obs : observations)
    {
        sigTotal.locus += obs.observationCount;
        unsigned totalAltObservations(0);
        for (unsigned altIndex(altBeginIndex); altIndex<altEndIndex; ++altIndex)
        {
            totalAltObservations += obs.altObservations[altIndex];
        }
        const unsigned total(obs.refObservations+totalAltObservations);
        const double altFrac(static_cast<double>(totalAltObservations)/total);
        if ( (total<minDepth) || (altFrac > maxAltFrac)) continue;

        sigTotal.ref += (obs.refObservations*obs.observationCount);
        sigTotal.alt += (totalAltObservations*obs.observationCount);
        sigTotal.noiseLocus += obs.observationCount;
    }

    {
        static const std::string sep(", ");
        const double total(sigTotal.ref+sigTotal.alt);

        static const double alpha(0.05);
        const double upper(boost::math::binomial_distribution<double>::find_upper_bound_on_p(total, sigTotal.alt, alpha));

        os << std::setprecision(10);
        os << context << "_" << extendedContextTag << sep
           << data.excludedRegionSkipped << sep
           << (sigTotal.locus + data.depthSkipped) << sep
           << sigTotal.locus << sep
           << sigTotal.noiseLocus << sep
           << sigTotal.ref << sep
           << sigTotal.alt << sep
           << safeFrac(sigTotal.alt,total) << sep
           << upper << "\n";
        ;
    }
}

}


void
indelModel1(
    const SequenceErrorCounts& counts)
{
    static const double maxAltFrac(0.05);
    static const unsigned minDepth(25);

    std::ostream& ros(std::cout);

    ros << "context, excludedLoci, nonExcludedLoci, usedLoci, noiseLoci, refReads, altReads, rate, rate_95%_upper_bound\n";

    std::vector<ExportedIndelObservations> observations;
    for (const auto& contextInfo : counts.getIndelCounts())
    {
        const auto& context(contextInfo.first);
        const auto& data(contextInfo.second);

        data.exportObservations(observations);

        if (observations.empty()) continue;

        reportExtendedContext(maxAltFrac, minDepth, context, observations, data,
                              INDEL_SIGNAL_TYPE::INSERT_1, INDEL_SIGNAL_TYPE::DELETE_1, "I", ros);
        reportExtendedContext(maxAltFrac, minDepth, context, observations, data,
                              INDEL_SIGNAL_TYPE::DELETE_1, INDEL_SIGNAL_TYPE::SIZE, "D", ros);
    }
}

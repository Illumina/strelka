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

#include "snvModel1.hh"

#include "blt_util/math_util.hh"
#include "blt_util/qscore.hh"

#include "boost/math/distributions/binomial.hpp"

#include <iomanip>
#include <iostream>
#include <numeric>
#include <set>



namespace
{

struct ContextTotal
{
    uint64_t locus = 0;
    uint64_t noiseLocus = 0;
};

}



template <typename K, typename V>
static
void
mergeMapKeys(
    const std::map<K,V>& m1,
    std::map<K,V>& m2,
    const unsigned scale = 1)
{
    for (const auto& mv1 : m1)
    {
        const auto iter(m2.find(mv1.first));
        if (iter == m2.end())
        {
            m2.insert(mv1);
        }
        else
        {
            iter->second += (mv1.second * scale);
        }
    }
}



/// print results for each context
///
static
void
printExtendedContext(
    const BaseErrorContext& context,
    const BaseErrorData& data,
    const ContextTotal& cTotal,
    const double refReadTotal,
    const double noiseLocusRefReadTotal,
    const StrandBaseCounts::qual_count_t& noiseLocusAltQualCounts,
    std::ostream& os)
{
    static const std::string sep(", ");
    const uint64_t emptySkippedLociCount(data.emptySkipped);
    const uint64_t depthSkippedLociCount(data.depthSkipped);
    const uint64_t noisyContextSkippedLociCount(data.noiseSkipped);

    const uint64_t totalNonExcludedLoci(cTotal.locus+emptySkippedLociCount+depthSkippedLociCount+noisyContextSkippedLociCount);

    const double refReadNoiseScaleFactor(noiseLocusRefReadTotal/refReadTotal);

    // get the set of qvals considered for this context:
    const auto& refQuals(data.error.getRefQuals());
    std::set<uint16_t> quals;
    for (const auto& value : data.error.getRefQuals())
    {
        quals.insert(value.first);
    }
    // TODO add quals form alt here:


    for (const auto qual : quals)
    {
        const double expect(qphred_to_error_prob(qual));
        const double refQualCount(refQuals.find(qual)->second);
        const double noiseLocusRefQualCount(refReadNoiseScaleFactor*refQualCount);
        const double noiseLocusAltQualCount(noiseLocusAltQualCounts.find(qual)->second);

        const double noiseLocusTotalQualCount(noiseLocusRefQualCount+noiseLocusAltQualCount);
        const double rate(safeFrac(noiseLocusAltQualCount,noiseLocusTotalQualCount));

        static const double alpha(0.05);
        const double upper(boost::math::binomial_distribution<double>::find_upper_bound_on_p(noiseLocusTotalQualCount, noiseLocusAltQualCount, alpha));

        os << context
           << sep << "Q" << qual
           << sep << data.excludedRegionSkipped
           << sep << totalNonExcludedLoci
           << sep << cTotal.locus
           << sep << cTotal.noiseLocus
           << sep << noiseLocusRefQualCount
           << sep << noiseLocusAltQualCount
           << sep << rate
           << sep << upper
           << sep << expect
           << "\n";
    }
}



static
void
reportExtendedContext(
    const double maxAltFrac,
    const unsigned minDepth,
    const BaseErrorContext& context,
    const BaseErrorData& data,
    std::ostream& os)
{
    //
    // Iterate through observations and categorize each as being possibly a variant
    // or an observation type we can use to estimate noise rates. While doing so,
    // compute aggregate statistics for context:
    // - locus count
    // - noise locus count
    // - total ref reads at noise loci
    // - total alt reads at noise loci
    //
    double refReadTotal(0);
    double noiseLocusRefReadTotal(0);
    StrandBaseCounts::qual_count_t noiseLocusAltQualCounts;

    ContextTotal cTotal;
    for (const auto& value : data.error)
    {
        const auto& key(value.first);
        const unsigned obsCount(value.second);

        const StrandBaseCounts& s0(key.getStrand0Counts());
        const StrandBaseCounts& s1(key.getStrand1Counts());

        const unsigned refQualTotal(s0.refCount+s1.refCount);

        cTotal.locus += obsCount;
        refReadTotal += refQualTotal*obsCount;

        unsigned altQualTotal(0);
        auto altSum = [&](const StrandBaseCounts& sb)
        {
            for (const auto& val : sb.alt)
            {
                altQualTotal += val.second;
            }
        };

        altSum(s0);
        altSum(s1);

        const unsigned depth(refQualTotal+altQualTotal);
        const double frac(altQualTotal/static_cast<double>(depth));
        if ((depth<minDepth) || (frac>maxAltFrac)) continue;

        cTotal.noiseLocus += obsCount;
        noiseLocusRefReadTotal += refQualTotal*obsCount;
        mergeMapKeys(s0.alt,noiseLocusAltQualCounts,obsCount);
        mergeMapKeys(s1.alt,noiseLocusAltQualCounts,obsCount);
    }

    printExtendedContext(context, data, cTotal, refReadTotal,
                         noiseLocusRefReadTotal, noiseLocusAltQualCounts, os);
}


void
snvModel1(
    const SequenceErrorCounts& counts)
{
    static const double maxAltFrac(0.03);
    static const unsigned minDepth(30);
    std::ostream& ros(std::cout);

    ros << "context, qual, excludedLoci, nonExcludedLoci, usedLoci, noiseLoci, refReads, altReads, rate, rate_95%_upper_bound, expect\n";

    for (const auto& contextInfo : counts.getBaseCounts())
    {
        reportExtendedContext(maxAltFrac,minDepth,contextInfo.first,contextInfo.second,ros);
    }
}

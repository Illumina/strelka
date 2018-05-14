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

#include "IndelCounts.hh"
#include "errorAnalysisUtils.hh"
#include "blt_util/IntegerLogCompressor.hh"
#include "blt_util/math_util.hh"

#include <cassert>
#include <cmath>

#include <iostream>


namespace IndelCounts
{

std::ostream&
operator<<(
    std::ostream& os,
    const Context& context)
{
    os << context.getRepeatPatternSize()
       << "x"
       << context.getRepeatCount();
    return os;
}

std::ostream&
operator<<(
    std::ostream& os,
    const SingleSampleNonVariantContextObservationPattern& obs)
{
    os << "Background observation\n";
    os << "refCount:" << obs.depth << std::endl;
    os << "variantStatus: " << GENOTYPE_STATUS::label(obs.backgroundStatus) << std::endl;
    return os;
}

std::ostream&
operator<<(
    std::ostream& os,
    const SingleSampleCandidateVariantContextObservationPattern& obs)
{
    os << "Indel observation\n";
    os << "refCount:" << obs.refCount << std::endl;
    os << "altCounts:\n";
    for (unsigned i(0); i < INDEL_SIGNAL_TYPE::SIZE; ++i)
    {
        os << "\t" << INDEL_SIGNAL_TYPE::label(i) << ":";
        os << obs.signalCounts[i] << std::endl;
    }
    os << "variantStatus: " << GENOTYPE_STATUS::label(obs.variantStatus) << std::endl;

    return os;
}

}

using namespace IndelCounts;



void
SingleSampleNonVariantContextData::
addSingleSampleNonVariantInstanceObservation(
    const SingleSampleNonVariantContextObservationPattern& obs)
{
    // one dimensional key doesn't need to be compressed as much (compared to the compression applied to the
    // more complex key used for variant counts)
    static const unsigned depthBitCount(6);

    SingleSampleNonVariantContextObservationPattern compObs = obs;
    compObs.depth = compressInt(compObs.depth,depthBitCount);

    iterateMapValue(data,compObs);
}



void
SingleSampleNonVariantContextData::
merge(const SingleSampleNonVariantContextData& in)
{
    mergeMapKeys(in.data,data);
}



void
SingleSampleNonVariantContextData::
dump(
    std::ostream& os) const
{
    uint64_t totalObservations(0.);
    uint64_t totalUnknownObservations(0);
    double totalDepth(0.);
    for (const auto& value : data)
    {
        totalObservations += value.second;
        totalDepth += (value.second*value.first.depth);
        if (value.first.backgroundStatus == GENOTYPE_STATUS::UNKNOWN)
        {
            totalUnknownObservations += value.second;
        }
    }

    const unsigned keyCount(data.size());

    static const std::string tag("background");
    os << tag << "KeyCount: " << keyCount << "\n";
    os << tag << "TotalObservations: " <<  totalObservations << "\n";
    os << tag << "TotalUnknownObservations: " << totalUnknownObservations << "\n";
    os << tag << "MeanKeyOccupancy: " <<  safeFrac(totalObservations,keyCount) << "\n";
    os << tag << "MeanDepth: " << safeFrac(totalDepth,totalObservations) << "\n";
}



void
SingleSampleCandidateVariantContextData::
addSingleSampleCandidateContextInstanceObservation(
    const SingleSampleCandidateVariantContextObservationPattern& obs)
{
    static const unsigned bitCount(5);

    SingleSampleCandidateVariantContextObservationPattern compObs = obs;
    compObs.refCount    = compressInt(compObs.refCount,bitCount);
    for (auto& signalCount : compObs.signalCounts)
    {
        signalCount = compressInt(signalCount,bitCount);
    }
    iterateMapValue(data,compObs);
}



void
SingleSampleCandidateVariantContextData::
merge(const SingleSampleCandidateVariantContextData& in)
{
    mergeMapKeys(in.data,data);
}



void
SingleSampleCandidateVariantContextData::
dump(
    std::ostream& os) const
{
    uint64_t totalObservations(0.);
    uint64_t totalUnknownObservations(0.);
    double totalRef(0.);
    double totalSignal(0.);

    unsigned noiseKeyCount(0);
    uint64_t noiseUnknownObservations(0.);
    uint64_t noiseObservations(0.);
    double noiseRef(0.);
    double noiseSignal(0.);

    for (const auto& value : data)
    {
        const auto& key(value.first);
        const auto& obsCount(value.second);
        totalObservations += obsCount;
        totalRef += (obsCount*key.refCount);
        totalSignal += (obsCount*key.totalSignalCount());
        if (key.variantStatus == GENOTYPE_STATUS::UNKNOWN)
        {
            totalUnknownObservations += obsCount;
        }

        const unsigned total(key.totalCount());
        const double frac(static_cast<double>(key.totalSignalCount())/total);
        if ((total >= 25) && (frac <= 0.05))
        {
            noiseKeyCount++;
            noiseObservations += obsCount;
            noiseRef += (obsCount*key.refCount);
            noiseSignal += (obsCount*key.totalSignalCount());
            if (key.variantStatus == GENOTYPE_STATUS::UNKNOWN)
            {
                noiseUnknownObservations += obsCount;
            }
        }
    }

    const unsigned keyCount(data.size());

    static const std::string tag("error");
    os << tag << "KeyCount: " << keyCount << "\n";
    os << tag << "TotalObservations: " <<  totalObservations << "\n";
    os << tag << "UnknownObservations: " << totalUnknownObservations << "\n";
    os << tag << "MeanKeyOccupancy: " <<  safeFrac(totalObservations,keyCount) << "\n";
    os << tag << "MeanRef: " << safeFrac(totalRef,totalObservations) << "\n";
    os << tag << "MeanSignal: " << safeFrac(totalSignal,totalObservations) << "\n";
    os << tag << "NoiseObservations: " <<  noiseObservations << "\n";
    os << tag << "NoiseUnknownObservations: " << noiseUnknownObservations << "\n";
    os << tag << "NoiseKeyCount: " << noiseKeyCount << "\n";
    os << tag << "MeanNoiseKeyOccupancy: " <<  safeFrac(noiseObservations,noiseKeyCount) << "\n";
    os << tag << "MeanNoiseRef: " << safeFrac(noiseRef,noiseObservations) << "\n";
    os << tag << "MeanNoiseSignal: " << safeFrac(noiseSignal,noiseObservations) << "\n";
}


namespace IndelCounts
{

std::ostream&
operator<<(
    std::ostream& os,
    const SingleSampleContextObservationInfoExportFormat& obs)
{
    os << "variantStatus" << GENOTYPE_STATUS::label(obs.variantStatus);
    os << "repeatCount: " << obs.contextInstanceCount;
    os << "\trefCounts: " << obs.refObservations;
    os << "\taltObservations: ";
    unsigned alt_index(0);
    for (const auto& alt: obs.altObservations)
    {
        os << INDEL_SIGNAL_TYPE::label(alt_index) << ":";
        os << alt << ",";
        ++alt_index;
    }

    return os;
}

}


void
ContextData::
exportData(
    SingleSampleContextDataExportFormat& exportedContextData) const
{
    exportedContextData.data.clear();

    // if this is true, we observed no errors:
    if (depthSupport.depth <= 0.) return;

    const double supportFraction(safeFrac(depthSupport.supportCount, depthSupport.depth));

    SingleSampleContextObservationInfoExportFormat contextObservationInfo;
    std::fill(contextObservationInfo.altObservations.begin(), contextObservationInfo.altObservations.end(), 0);

    // convert background observations:
    for (const auto& value : nonVariantContextCounts)
    {
        contextObservationInfo.contextInstanceCount = value.second;
        contextObservationInfo.refObservations = (value.first.depth*supportFraction);
        contextObservationInfo.variantStatus = value.first.backgroundStatus;

        exportedContextData.data.push_back(contextObservationInfo);
    }

    // convert error observations:
    for (const auto& value : candidateVariantContextCounts)
    {
        contextObservationInfo.contextInstanceCount = value.second;
        contextObservationInfo.refObservations = value.first.refCount;
        contextObservationInfo.altObservations = value.first.signalCounts;
        contextObservationInfo.variantStatus   = value.first.variantStatus;

        exportedContextData.data.push_back(contextObservationInfo);
    }
}

void
ContextData::
merge(
    const ContextData& in)
{
    nonVariantContextCounts.merge(in.nonVariantContextCounts);
    candidateVariantContextCounts.merge(in.candidateVariantContextCounts);
    depthSupport.merge(in.depthSupport);
    excludedRegionSkipped += in.excludedRegionSkipped;
    depthSkipped += in.depthSkipped;
}



void
ContextData::
dump(
    std::ostream& os) const
{
    os << "depthSupportRatio: " << safeFrac(depthSupport.supportCount, depthSupport.depth)
       << " (" << depthSupport.supportCount << "/" << depthSupport.depth << ")\n";
    os << "excludedRegionSkippedCount: " << excludedRegionSkipped << "\n";
    os << "depthSkippedCount: " << depthSkipped << "\n";

    nonVariantContextCounts.dump(os);
    candidateVariantContextCounts.dump(os);
}



void
Dataset::
addCandidateVariantContextInstanceObservation(
    const Context& context,
    const SingleSampleCandidateVariantContextObservationPattern& candidateVariantContextObservation,
    const unsigned depth)
{
    const auto iter(getContextIterator(context));
    iter->second.candidateVariantContextCounts.addSingleSampleCandidateContextInstanceObservation(candidateVariantContextObservation);
    iter->second.depthSupport.merge(IndelDepthSupportTotal(depth,candidateVariantContextObservation.totalCount()));
}



void
Dataset::
addNonVariantContextInstanceObservation(
    const Context& context,
    const SingleSampleNonVariantContextObservationPattern& nonVariantContextObservation)
{
    const auto iter(getContextIterator(context));
    iter->second.nonVariantContextCounts.addSingleSampleNonVariantInstanceObservation(nonVariantContextObservation);
}



void
Dataset::
addExcludedRegionSkip(
    const Context& context)
{
    const auto iter(getContextIterator(context));
    iter->second.excludedRegionSkipped++;
}



void
Dataset::
addDepthSkip(
    const Context& context)
{
    const auto iter(getContextIterator(context));
    iter->second.depthSkipped++;
}



void
Dataset::
merge(
    const Dataset& in)
{
    for (const auto& idata : in._data)
    {
        const auto iter(_data.find(idata.first));

        if (iter == _data.end())
        {
            _data.insert(idata);
        }
        else
        {
            iter->second.merge(idata.second);
        }
    }
}



void
Dataset::
dump(
    std::ostream& os) const
{
    os << "IndelCounts DUMP_ON\n";
    os << "Total Indel Contexts: " << _data.size() << "\n";
    for (const auto& value : _data)
    {
        os << "Indel Context: Repeat pattern size " << value.first.getRepeatPatternSize() << ", Repeat count " << value.first.getRepeatCount() << "\n";
        value.second.dump(os);
    }
    os << "IndelCounts DUMP_OFF\n";
}



Dataset::data_t::iterator
Dataset::
getContextIterator(
    const Context& context)
{
    const auto insert = _data.insert(std::make_pair(context,data_t::mapped_type()));
    return insert.first;
}

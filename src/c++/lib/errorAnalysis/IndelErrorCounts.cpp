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

#include "IndelErrorCounts.hh"
#include "blt_util/IntegerLogCompressor.hh"
#include "blt_util/math_util.hh"

#include <cassert>
#include <cmath>

#include <iostream>



/// assuming V is an integer count type, iterate or initialize
/// new key
template <typename K, typename V>
static
void
iterMap(
    std::map<K,V>& m,
    const K& key)
{
    const auto iter(m.find(key));
    if (iter == m.end())
    {
        m[key] = 1;
    }
    else
    {
        iter->second += 1;
    }
}


template <typename K, typename V>
static
void
mergeMapKeys(
    const std::map<K,V>& m1,
    std::map<K,V>& m2)
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
            iter->second += mv1.second;
        }
    }
}



std::ostream&
operator<<(
    std::ostream& os,
    const IndelErrorContext& context)
{
    os << context.getRepeatPatternSize()
       << "x"
       << context.getRepeatCount();
    return os;
}

std::ostream&
operator<<(
    std::ostream& os,
    const IndelBackgroundObservation& obs)
{
    os << "Background observation\n";
    os << "refCount:" << obs.depth << std::endl;
    os << "variantStatus: " << GENOTYPE_STATUS::label(obs.backgroundStatus) << std::endl;
    return os;
}

std::ostream&
operator<<(
    std::ostream& os,
    const IndelErrorContextObservation& obs)
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

void
IndelBackgroundObservationData::
addObservation(
    const IndelBackgroundObservation& obs)
{
    // 1D key doesn't need to be compressed as much:
    static const unsigned depthBitCount(6);

    IndelBackgroundObservation compObs = obs;
    compObs.depth = compressInt(compObs.depth,depthBitCount);

    iterMap(data,compObs);
}



void
IndelBackgroundObservationData::
merge(const IndelBackgroundObservationData& in)
{
    mergeMapKeys(in.data,data);
}



void
IndelBackgroundObservationData::
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
IndelErrorContextObservationData::
addObservation(
    const IndelErrorContextObservation& obs)
{
    static const unsigned bitCount(5);

    IndelErrorContextObservation compObs = obs;
    compObs.refCount    = compressInt(compObs.refCount,bitCount);
    for (auto& signalCount : compObs.signalCounts)
    {
        signalCount = compressInt(signalCount,bitCount);
    }
    iterMap(data,compObs);
}



void
IndelErrorContextObservationData::
merge(const IndelErrorContextObservationData& in)
{
    mergeMapKeys(in.data,data);
}



void
IndelErrorContextObservationData::
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



std::ostream&
operator<<(
    std::ostream& os,
    const ExportedIndelObservations& obs)
{
    os << "variantStatus" << GENOTYPE_STATUS::label(obs.variantStatus);
    os << "repeatCount: " << obs.observationCount;
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



void
IndelErrorData::
exportObservations(
    std::vector<ExportedIndelObservations>& counts) const
{
    counts.clear();

    // if this is true, we observed no errors:
    if (depthSupport.depth <= 0.) return;

    const double supportFraction(safeFrac(depthSupport.supportCount, depthSupport.depth));

    ExportedIndelObservations obs;
    std::fill(obs.altObservations.begin(), obs.altObservations.end(), 0);

    // convert background observations:
    for (const auto& value : background)
    {
        obs.observationCount = value.second;
        obs.refObservations = (value.first.depth*supportFraction);
        obs.variantStatus = value.first.backgroundStatus;

        counts.push_back(obs);
    }

    // convert error observations:
    for (const auto& value : error)
    {
        obs.observationCount = value.second;
        obs.refObservations = value.first.refCount;
        obs.altObservations = value.first.signalCounts;
        obs.variantStatus   = value.first.variantStatus;

        counts.push_back(obs);
    }
}

void
IndelErrorData::
merge(
    const IndelErrorData& in)
{
    background.merge(in.background);
    error.merge(in.error);
    depthSupport.merge(in.depthSupport);
    excludedRegionSkipped += in.excludedRegionSkipped;
    depthSkipped += in.depthSkipped;
}



void
IndelErrorData::
dump(
    std::ostream& os) const
{
    os << "depthSupportRatio: " << safeFrac(depthSupport.supportCount, depthSupport.depth)
       << " (" << depthSupport.supportCount << "/" << depthSupport.depth << ")\n";
    os << "excludedRegionSkippedCount: " << excludedRegionSkipped << "\n";
    os << "depthSkippedCount: " << depthSkipped << "\n";

    background.dump(os);
    error.dump(os);
}



void
IndelErrorCounts::
addError(
    const IndelErrorContext& context,
    const IndelErrorContextObservation& errorObservation,
    const unsigned depth)
{
    const auto iter(getContextIterator(context));
    iter->second.error.addObservation(errorObservation);
    iter->second.depthSupport.merge(IndelDepthSupportTotal(depth,errorObservation.totalCount()));
}



void
IndelErrorCounts::
addBackground(
    const IndelErrorContext& context,
    const IndelBackgroundObservation& backgroundObservation)
{
    const auto iter(getContextIterator(context));
    iter->second.background.addObservation(backgroundObservation);
}



void
IndelErrorCounts::
addExcludedRegionSkip(
    const IndelErrorContext& context)
{
    const auto iter(getContextIterator(context));
    iter->second.excludedRegionSkipped++;
}



void
IndelErrorCounts::
addDepthSkip(
    const IndelErrorContext& context)
{
    const auto iter(getContextIterator(context));
    iter->second.depthSkipped++;
}



void
IndelErrorCounts::
merge(
    const IndelErrorCounts& in)
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
IndelErrorCounts::
dump(
    std::ostream& os) const
{
    os << "IndelErrorCounts DUMP_ON\n";
    os << "Total Indel Contexts: " << _data.size() << "\n";
    for (const auto& value : _data)
    {
        os << "Indel Context: Repeat pattern size " << value.first.getRepeatPatternSize() << ", Repeat count " << value.first.getRepeatCount() << "\n";
        value.second.dump(os);
    }
    os << "IndelErrorCounts DUMP_OFF\n";
}



IndelErrorCounts::data_t::iterator
IndelErrorCounts::
getContextIterator(
    const IndelErrorContext& context)
{
    const auto insert = _data.insert(std::make_pair(context,data_t::mapped_type()));
    return insert.first;
}

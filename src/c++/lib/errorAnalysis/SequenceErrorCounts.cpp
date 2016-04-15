// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2016 Illumina, Inc.
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

#include "SequenceErrorCounts.hh"
#include "blt_util/IntegerLogCompressor.hh"
#include "blt_util/math_util.hh"

#include "boost/archive/binary_iarchive.hpp"
#include "boost/archive/binary_oarchive.hpp"
#include "boost/serialization/array.hpp"
#include "boost/serialization/map.hpp"

#include <cassert>
#include <cmath>

#include <fstream>
#include <iostream>


/// assuming V is an integer count type, iterate or initialize
/// new key
template <typename K, typename V>
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
    const SequenceErrorContext& context)
{
    os << context.repeatCount;
    return os;
}



void
SequenceBackgroundObservationData::
addObservation(
    const SequenceBackgroundObservation& obs)
{
    // 1D key doesn't need to be compressed as much:
    static const unsigned depthBitCount(6);

    SequenceBackgroundObservation compObs = obs;
    compObs.depth = compressInt(compObs.depth,depthBitCount);

    iterMap(data,compObs);
}



void
SequenceBackgroundObservationData::
merge(const SequenceBackgroundObservationData& in)
{
    mergeMapKeys(in.data,data);
}



void
SequenceBackgroundObservationData::
dump(
    std::ostream& os) const
{
    uint64_t totalObservations(0.);
    double totaDepth(0.);
    for (const auto& value : data)
    {
        totalObservations += value.second;
        totaDepth += (value.second*value.first.depth);
    }

    const unsigned keyCount(data.size());

    static const std::string tag("background");
    os << tag << "KeyCount: " << keyCount << "\n";
    os << tag << "TotalObservations: " <<  totalObservations << "\n";
    os << tag << "MeanKeyOccupancy: " <<  safeFrac(totalObservations,keyCount) << "\n";
    os << tag << "MeanDepth: " << safeFrac(totaDepth,totalObservations) << "\n";
}



void
SequenceErrorContextObservationData::
addObservation(
    const SequenceErrorContextObservation& obs)
{
    static const unsigned bitCount(5);

    SequenceErrorContextObservation compObs = obs;
    compObs.refCount    = compressInt(compObs.refCount,bitCount);
    for (auto& signalCount : compObs.signalCounts)
    {
        signalCount = compressInt(signalCount,bitCount);
    }
    iterMap(data,compObs);
}



void
SequenceErrorContextObservationData::
merge(const SequenceErrorContextObservationData& in)
{
    mergeMapKeys(in.data,data);
}



void
SequenceErrorContextObservationData::
dump(
    std::ostream& os) const
{
    uint64_t totalObservations(0.);
    double totalRef(0.);
    double totalSignal(0.);

    unsigned noiseKeyCount(0);
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

        const unsigned total(key.totalCount());
        const double frac(static_cast<double>(key.totalSignalCount())/total);
        if ((total >= 25) && (frac <= 0.05))
        {
            noiseKeyCount++;
            noiseObservations += obsCount;
            noiseRef += (obsCount*key.refCount);
            noiseSignal += (obsCount*key.totalSignalCount());
        }
    }

    const unsigned keyCount(data.size());

    static const std::string tag("error");
    os << tag << "KeyCount: " << keyCount << "\n";
    os << tag << "TotalObservations: " <<  totalObservations << "\n";
    os << tag << "MeanKeyOccupancy: " <<  safeFrac(totalObservations,keyCount) << "\n";
    os << tag << "MeanRef: " << safeFrac(totalRef,totalObservations) << "\n";
    os << tag << "MeanSignal: " << safeFrac(totalSignal,totalObservations) << "\n";
    os << tag << "NoiseObservations: " <<  noiseObservations << "\n";
    os << tag << "NoiseKeyCount: " << noiseKeyCount << "\n";
    os << tag << "MeanNoiseKeyOccupancy: " <<  safeFrac(noiseObservations,noiseKeyCount) << "\n";
    os << tag << "MeanNoiseRef: " << safeFrac(noiseRef,noiseObservations) << "\n";
    os << tag << "MeanNoiseSignal: " << safeFrac(noiseSignal,noiseObservations) << "\n";
}

std::ostream&
operator<<(
    std::ostream& os,
    const ExportedObservations& obs)
{
    os << "repeatCount: " << obs.repeatCount;
    os << "\trefCounts: " << obs.refObservations;
    os << "\taltObservations: ";
    unsigned alt_index(0);
    for (const auto& alt: obs.altObservations)
    {
        os << SIGNAL_TYPE::label(static_cast<SIGNAL_TYPE::index_t>(alt_index)) << ":";
        os << alt << ",";
        ++alt_index;
    }

    return os;
}


void
SequenceErrorData::
exportObservations(
    std::vector<ExportedObservations>& counts) const
{
    counts.clear();

    // if this is true, we observed no errors:
    if (depthSupport.depth <= 0.) return;

    const double supportFraction(safeFrac(depthSupport.supportCount, depthSupport.depth));

    ExportedObservations obs;

    // convert background observations:
    for (const auto& value : background)
    {
        const unsigned refCounts(std::round(value.first.depth*supportFraction));
        obs.repeatCount = value.second;
        obs.refObservations = refCounts;
        std::fill(obs.altObservations.begin(), obs.altObservations.end(), 0);

        counts.push_back(obs);
    }

    // convert error observations:
    for (const auto& value : error)
    {
        obs.repeatCount = value.second;
        obs.refObservations = value.first.refCount;
        obs.altObservations = value.first.signalCounts;

        counts.push_back(obs);
    }
}



void
SequenceErrorData::
merge(
    const SequenceErrorData& in)
{
    background.merge(in.background);
    error.merge(in.error);
    depthSupport.merge(in.depthSupport);
    skipped += in.skipped;
}



void
SequenceErrorData::
dump(
    std::ostream& os) const
{
    os << "depthSupportRatio: " << safeFrac(depthSupport.supportCount, depthSupport.depth) << "\n";
    os << "skippedCount: " << skipped << "\n";

    background.dump(os);
    error.dump(os);
}



void
SequenceErrorCounts::
addError(
    const SequenceErrorContext& context,
    const SequenceErrorContextObservation& errorObservation,
    const unsigned depth)
{
    const auto iter(getContextIterator(context));
    iter->second.error.addObservation(errorObservation);
    iter->second.depthSupport.merge(SequenceDepthSupportTotal(depth,errorObservation.totalCount()));
}



void
SequenceErrorCounts::
addBackground(
    const SequenceErrorContext& context,
    const SequenceBackgroundObservation& backgroundObservation)
{
    const auto iter(getContextIterator(context));
    iter->second.background.addObservation(backgroundObservation);
}



void
SequenceErrorCounts::
addDepthSkip(
    const SequenceErrorContext& context)
{
    const auto iter(getContextIterator(context));
    iter->second.skipped++;
}



void
SequenceErrorCounts::
merge(
    const SequenceErrorCounts& in)
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
SequenceErrorCounts::
save(
    const char* filename) const
{
    using namespace boost::archive;

    assert(nullptr != filename);
    std::ofstream ofs(filename, std::ios::binary);
    binary_oarchive oa(ofs);

    oa << _data;
}



void
SequenceErrorCounts::
load(
    const char* filename)
{
    using namespace boost::archive;

    clear();

    assert(nullptr != filename);
    std::ifstream ifs(filename, std::ios::binary);
    binary_iarchive ia(ifs);

    ia >> _data;
}



void
SequenceErrorCounts::
dump(
    std::ostream& os) const
{
    os << "SequenceErrorCounts DUMP_ON\n";
    os << "Total Contexts: " << _data.size() << "\n";
    for (const auto& value : _data)
    {
        os << "Context: " << value.first << "\n";
        value.second.dump(os);
    }
    os << "SequenceErrorCounts DUMP_OFF\n";
}


SequenceErrorCounts::data_t::iterator
SequenceErrorCounts::
getContextIterator(
    const SequenceErrorContext& context)
{
    const auto insert = _data.insert(std::make_pair(context,data_t::mapped_type()));
    return insert.first;
}

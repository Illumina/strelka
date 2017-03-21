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

#include "BaseErrorCounts.hh"
#include "blt_util/IntegerLogCompressor.hh"
#include "blt_util/math_util.hh"

#include <cassert>
#include <cmath>

#include <iostream>
#include <set>



template <typename K, typename V>
static
V
mapFindDefault(
    const std::map<K,V>& m,
    const K& key,
    const V& defaultVal)
{
    const auto iter(m.find(key));
    if (iter == m.end())
    {
        return defaultVal;
    }
    else
    {
        return iter->second;
    }
}



/// assuming V is an integer count type, iterate or initialize
/// new key
template <typename K, typename V>
static
void
iterMap(
    std::map<K,V>& m,
    const K& key,
    const unsigned iterValue = 1)
{
    const auto iter(m.find(key));
    if (iter == m.end())
    {
        m[key] = iterValue;
    }
    else
    {
        iter->second += iterValue;
    }
}



template <typename K, typename V1, typename V2>
static
void
mergeMapKeys(
    const std::map<K,V1>& m1,
    std::map<K,V2>& m2,
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



std::ostream&
operator<<(
    std::ostream& os,
    const BaseErrorContext& context)
{
    os << context.repeatCount;
    return os;
}



void
StrandBaseCounts::
compressCounts()
{
    static const unsigned bitCount(4);

    refCount = compressInt(refCount, bitCount);
    for (auto& val : alt)
    {
        val.second = compressInt(val.second, bitCount);
    }
}



std::ostream&
operator<<(std::ostream& os, const StrandBaseCounts& sbc)
{
    os << "REF:\t" << sbc.refCount;

    bool isFirst=true;
    os << "\tALT:\t";
    for (const auto& val : sbc.alt)
    {
        if (isFirst)
        {
            isFirst=false;
        }
        else
        {
            os << ',';
        }
        os << val.first << ':' << val.second;
    }
    return os;
}



void
BaseErrorContextInputObservation::
addRefCount(
    const bool isFwdStrand,
    const uint16_t qual)
{
    StrandBaseCounts::qual_count_t& target(isFwdStrand ? ref[0] : ref[1]);
    iterMap(target,qual);
}



void
BaseErrorContextInputObservation::
addAltCount(
    const bool isFwdStrand,
    const uint16_t qual)
{
    StrandBaseCounts::qual_count_t& target(isFwdStrand ? alt[0] : alt[1]);
    iterMap(target,qual);
}



void
BaseErrorContextObservation::
compressCounts()
{
    // if no alts exist, we can safely erase strand information by summing everything to strand1:
    if (strand0.alt.empty() && strand1.alt.empty())
    {
        strand0.refCount += strand1.refCount;
        strand1.refCount=0;
    }

    strand0.compressCounts();
    strand1.compressCounts();
    if (strand0 < strand1)
    {
        std::swap(strand0,strand1);
    }
}



void
BaseErrorContextObservationData::
addObservation(
    const BaseErrorContextInputObservation& obs)
{
    BaseErrorContextObservation compObs;
    for (unsigned strandId(0); strandId<2; ++strandId)
    {
        mergeMapKeys(obs.ref[strandId], refQuals);
        auto& target(strandId==0 ? compObs.strand0 : compObs.strand1);
        for (const auto& val : obs.ref[strandId])
        {
            target.refCount += val.second;
        }
        target.alt = obs.alt[strandId];
    }

    compObs.compressCounts();
    iterMap(data,compObs);
}



void
BaseErrorContextObservationData::
merge(const BaseErrorContextObservationData& in)
{
    mergeMapKeys(in.data,data);
    mergeMapKeys(in.refQuals,refQuals);
}



void
BaseErrorContextObservationData::
getExportData(BaseErrorContextObservationExportData& exportData) const
{
    exportData.clear();

    //
    // find and set the final qual level list:
    //
    std::set<uint16_t> quals;

    // add qual levels from ref:
    for (const auto& val : refQuals)
    {
        quals.insert(val.first);
    }
    // ...don't even bother adding qual levels from alt, if there's an alt only level -- skip it, or assert.

    // tmp data structure for this function only:
    std::map<uint16_t,unsigned> qualIndex;

    for (const uint16_t qual : quals)
    {
        qualIndex[qual] = exportData.qualLevels.size();
        exportData.qualLevels.push_back(qual);
    }

    // set refCount wrt the full qual list:
    for (const uint16_t qual : exportData.qualLevels)
    {
        uint64_t refCount(0);
        const auto iter(refQuals.find(qual));
        if (iter != refQuals.end())
        {
            refCount = iter->second;
        }
        exportData.refCount.push_back(refCount);
    }

    // convert observations to export observations:
    for (const auto& value : data)
    {
        const BaseErrorContextObservation& key(value.first);
        const auto& obsCount(value.second);

        BaseErrorContextObservationExportObservation exportKey;

        auto strand2Strand = [&](
                                 const StrandBaseCounts& si,
                                 BaseErrorContextObservationExportStrandObservation& se)
        {
            se.refCount = si.refCount;
            se.altCount.resize(quals.size(),0);
            for (const auto& altValue : si.alt)
            {
                const uint16_t qual(altValue.first);
                const unsigned count(altValue.second);
                se.altCount[qualIndex.find(qual)->second] = count;
            }
        };

        strand2Strand(key.getStrand0Counts(), exportKey.strand0);
        strand2Strand(key.getStrand1Counts(), exportKey.strand1);

        assert(exportData.observations.find(exportKey) == exportData.observations.end());
        exportData.observations[exportKey] = obsCount;
    }
}



void
BaseErrorContextObservationData::
dump(
    std::ostream& os) const
{
    const unsigned keyCount(data.size());

    static const std::string tag("base-error");
    os << tag << "KeyCount: " << keyCount << "\n";

    unsigned refOnlyKeyCount(0);
    unsigned altOnlyKeyCount(0);
    uint64_t totalObservations(0.);
    refQual_t totalRef(refQuals);
    refQual_t totalAlt;

    std::map<unsigned,uint64_t> totalByDepth;

    for (const auto& value : data)
    {
        const BaseErrorContextObservation& key(value.first);
        const auto& obsCount(value.second);
        totalObservations += obsCount;

        const auto& s0(key.getStrand0Counts());
        const auto& s1(key.getStrand1Counts());
        mergeMapKeys(s0.alt,totalAlt,obsCount);
        mergeMapKeys(s1.alt,totalAlt,obsCount);

        // update depth map:
        {
            unsigned depth(0);
            depth += (s0.refCount + s1.refCount);
            for (const auto& altVal : s0.alt)
            {
                depth += altVal.second;
            }
            for (const auto& altVal : s1.alt)
            {
                depth += altVal.second;
            }
            iterMap(totalByDepth,depth,obsCount);
        }

        if (key.getStrand0Counts().alt.empty() && key.getStrand1Counts().alt.empty())
        {
            refOnlyKeyCount++;
        }
        if (key.getStrand0Counts().refCount==0 && key.getStrand1Counts().refCount==0)
        {
            altOnlyKeyCount++;
        }
    }

    os << tag << "RefOnlyKeyCount: " << refOnlyKeyCount << "\n";
    os << tag << "AltOnlyKeyCount: " << altOnlyKeyCount << "\n";
    os << tag << "TotalObservations: " <<  totalObservations << "\n";
    os << tag << "MeanKeyOccupancy: " <<  safeFrac(totalObservations,keyCount) << "\n";

    // get union of qual values form ref/alt:
    std::set<uint16_t> quals;
    for (const auto& value : totalRef)
    {
        quals.insert(value.first);
    }
    for (const auto& value : totalAlt)
    {
        quals.insert(value.first);
    }

    os << tag << "Qval\tTotalRef\tTotalAlt\n";

    for (const auto qual : quals)
    {

        const auto refCount(mapFindDefault(totalRef,qual,refQual_t::mapped_type(0)));
        const auto altCount(mapFindDefault(totalAlt,qual,refQual_t::mapped_type(0)));
        os << tag << "Q" << qual << "\t" << refCount << "\t" << altCount << "\n";
    }

    for (const auto& value : totalByDepth)
    {
        const unsigned depth(value.first);
        const uint64_t observations(value.second);
        os << "DEPTH: " << depth << "\t" << observations << "\n";
    }

    // uncomment this section to fully enumerate the SNV counts data:
#if 0
    for (const auto& value : data)
    {
        const BaseErrorContextObservation& key(value.first);
        const auto& obsCount(value.second);
        os << "KEYcount: " << obsCount << "\n";
        os << "KEYS1:\t" << key.getStrand0Counts() << "\n";
        os << "KEYS2:\t" << key.getStrand1Counts() << "\n";
    }
#endif
}



void
BaseErrorData::
merge(
    const BaseErrorData& in)
{
    error.merge(in.error);
    excludedRegionSkipped += in.excludedRegionSkipped;
    depthSkipped += in.depthSkipped;
    emptySkipped += in.emptySkipped;
    noiseSkipped += in.noiseSkipped;
}



void
BaseErrorData::
dump(
    std::ostream& os) const
{
    os << "excludedRegionSkippedCount: " << excludedRegionSkipped << "\n";
    os << "depthSkippedCount: " << depthSkipped << "\n";
    os << "emptySkippedCount: " << emptySkipped << "\n";
    os << "noiseSkippedCount: " << noiseSkipped << "\n";

    error.dump(os);
}



void
BaseErrorCounts::
addSiteObservation(
    const BaseErrorContext& context,
    const BaseErrorContextInputObservation& siteObservation)
{
    const auto iter(getContextIterator(context));
    iter->second.error.addObservation(siteObservation);
}



void
BaseErrorCounts::
addExcludedRegionSkip(
    const BaseErrorContext& context)
{
    const auto iter(getContextIterator(context));
    iter->second.excludedRegionSkipped++;
}



void
BaseErrorCounts::
addDepthSkip(
    const BaseErrorContext& context)
{
    const auto iter(getContextIterator(context));
    iter->second.depthSkipped++;
}



void
BaseErrorCounts::
addEmptySkip(
    const BaseErrorContext& context)
{
    const auto iter(getContextIterator(context));
    iter->second.emptySkipped++;
}



void
BaseErrorCounts::
addNoiseSkip(
    const BaseErrorContext& context)
{
    const auto iter(getContextIterator(context));
    iter->second.noiseSkipped++;
}



void
BaseErrorCounts::
merge(
    const BaseErrorCounts& in)
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
BaseErrorCounts::
dump(
    std::ostream& os) const
{
    os << "BaseErrorCounts DUMP_ON\n";
    os << "Total Base Contexts: " << _data.size() << "\n";
    for (const auto& value : _data)
    {
        os << "Base Context: " << value.first << "\n";
        value.second.dump(os);
    }
    os << "BaseErrorCounts DUMP_OFF\n";
}



BaseErrorCounts::data_t::iterator
BaseErrorCounts::
getContextIterator(
    const BaseErrorContext& context)
{
    const auto insert = _data.insert(std::make_pair(context,data_t::mapped_type()));
    return insert.first;
}

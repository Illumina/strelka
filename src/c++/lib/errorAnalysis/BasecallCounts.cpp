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


#include "BasecallCounts.hh"
#include "errorAnalysisUtils.hh"
#include "blt_util/IntegerLogCompressor.hh"
#include "blt_util/math_util.hh"

#include <cassert>
#include <cmath>

#include <iostream>
#include <set>


using namespace BasecallCounts;



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


namespace BasecallCounts
{

std::ostream&
operator<<(
    std::ostream& os,
    const Context& context)
{
    os << context.repeatCount;
    return os;
}

}


void
SingleSampleSingleStrandContextObservationPattern::
compressCounts()
{
    static const unsigned bitCount(4);

    refAlleleCount = compressInt(refAlleleCount, bitCount);
    for (auto& keyVal : altAlleleCount)
    {
        keyVal.second = compressInt(keyVal.second, bitCount);
    }
}


namespace BasecallCounts
{

std::ostream&
operator<<(
    std::ostream& os,
    const SingleSampleSingleStrandContextObservationPattern& pattern)
{
    os << "REF:\t" << pattern.refAlleleCount;

    bool isFirst = true;
    os << "\tALT:\t";
    for (const auto& val : pattern.altAlleleCount)
    {
        if (isFirst)
        {
            isFirst = false;
        }
        else
        {
            os << ',';
        }
        os << val.first << ':' << val.second;
    }
    return os;
}

}


void
ContextInstanceObservation::
addRefCount(
    const bool isFwdStrand,
    const uint16_t basecallErrorPhredProb)
{
    SingleSampleSingleStrandContextObservationPattern::qual_count_t& target(isFwdStrand ? ref[0] : ref[1]);
    iterateMapValue(target, basecallErrorPhredProb);
}



void
ContextInstanceObservation::
addAltCount(
    const bool isFwdStrand,
    const uint16_t basecallErrorPhredProb)
{
    SingleSampleSingleStrandContextObservationPattern::qual_count_t& target(isFwdStrand ? alt[0] : alt[1]);
    iterateMapValue(target, basecallErrorPhredProb);
}



void
SingleSampleContextObservationPattern::
compressCounts()
{
    // if no alts exist, we can safely erase strand information by summing everything to strand0
    // TODO: elaborate on why this is 'safe'
    if (strand0.altAlleleCount.empty() && strand1.altAlleleCount.empty())
    {
        strand0.refAlleleCount += strand1.refAlleleCount;
        strand1.refAlleleCount = 0;
    }

    strand0.compressCounts();
    strand1.compressCounts();
    if (strand0 < strand1)
    {
        std::swap(strand0,strand1);
    }
}



void
SingleSampleContextData::
addSingleSampleContextInstanceObservation(
    const ContextInstanceObservation& obs)
{
    static const unsigned strandCount(2);

    SingleSampleContextObservationPattern compressedObs;
    for (unsigned strandIndex(0); strandIndex<strandCount; ++strandIndex)
    {
        mergeMapKeys(obs.ref[strandIndex], refAlleleBasecallErrorPhredProbs);
        auto& targetStrand(strandIndex==0 ? compressedObs.strand0 : compressedObs.strand1);
        for (const auto& val : obs.ref[strandIndex])
        {
            // quality of reference allele observations are removed in the compressed observations
            targetStrand.refAlleleCount += val.second;
        }
        targetStrand.altAlleleCount = obs.alt[strandIndex];
    }

    compressedObs.compressCounts();
    iterateMapValue(data, compressedObs);
}



void
SingleSampleContextData::
merge(const SingleSampleContextData& in)
{
    mergeMapKeys(in.data,data);
    mergeMapKeys(in.refAlleleBasecallErrorPhredProbs,refAlleleBasecallErrorPhredProbs);
}



void
SingleSampleContextData::
exportData(SingleSampleContextDataExportFormat& exportedData) const
{
    exportedData.clear();

    //
    // find and the final set of basecall error levels:
    //
    std::set<uint16_t> basecallErrorPhredProbs;

    // add basecall error levels from reference allele observations:
    for (const auto& keyVal : refAlleleBasecallErrorPhredProbs)
    {
        basecallErrorPhredProbs.insert(keyVal.first);
    }
    // ...don't even bother adding basecall error levels from the alt alleles, it should be incredibly rare that a
    // level would be exclusive to the alt alleles. If this does occur, skip or assert.

    // tmp data structure used in this function only to map quality values to the corresponding quality array index
    std::map<uint16_t,unsigned> qualIndex;

    for (const uint16_t qual : basecallErrorPhredProbs)
    {
        qualIndex[qual] = exportedData.altAlleleBasecallErrorPhredProbLevels.size();
        exportedData.altAlleleBasecallErrorPhredProbLevels.push_back(qual);
    }

    // set refCount wrt the full qual list:
    for (const uint16_t qual : exportedData.altAlleleBasecallErrorPhredProbLevels)
    {
        uint64_t refCount(0);
        const auto iter(refAlleleBasecallErrorPhredProbs.find(qual));
        if (iter != refAlleleBasecallErrorPhredProbs.end())
        {
            refCount = iter->second;
        }
        exportedData.refCount.push_back(refCount);
    }

    // convert observations to export observations:
    for (const auto& keyVal : data)
    {
        const SingleSampleContextObservationPattern& observationPattern(keyVal.first);
        const auto& instanceCount(keyVal.second);

        SingleSampleContextObservationPatternExportFormat exportObservationPattern;

        /// Convert single strand observation counts from compressed storage format to exported format
        /// intended for use by an inference method
        auto strand2ExportStrand = [&](
                                       const SingleSampleSingleStrandContextObservationPattern& si,
                                       SingleSampleSingleStrandContextObservationPatternExportFormat& se)
        {
            se.refAlleleCount = si.refAlleleCount;
            se.altAlleleCount.resize(basecallErrorPhredProbs.size(),0);
            for (const auto& qualCount : si.altAlleleCount)
            {
                se.altAlleleCount[qualIndex.find(qualCount.first)->second] = qualCount.second;
            }
        };

        strand2ExportStrand(observationPattern.getStrand0Counts(), exportObservationPattern.strand0);
        strand2ExportStrand(observationPattern.getStrand1Counts(), exportObservationPattern.strand1);

        assert(exportedData.observations.find(exportObservationPattern) == exportedData.observations.end());
        exportedData.observations[exportObservationPattern] = instanceCount;
    }
}



void
SingleSampleContextData::
dump(
    std::ostream& os) const
{
    const unsigned keyCount(data.size());

    static const std::string tag("base-error");
    os << tag << "KeyCount: " << keyCount << "\n";

    unsigned refOnlyKeyCount(0);
    unsigned altOnlyKeyCount(0);
    uint64_t totalObservations(0.);
    refQual_t totalRef(refAlleleBasecallErrorPhredProbs);
    refQual_t totalAlt;

    std::map<unsigned,uint64_t> totalByDepth;

    for (const auto& value : data)
    {
        const SingleSampleContextObservationPattern& key(value.first);
        const auto& obsCount(value.second);
        totalObservations += obsCount;

        const auto& s0(key.getStrand0Counts());
        const auto& s1(key.getStrand1Counts());
        mergeMapKeys(s0.altAlleleCount,totalAlt,obsCount);
        mergeMapKeys(s1.altAlleleCount,totalAlt,obsCount);

        // update depth map:
        {
            unsigned depth(0);
            depth += (s0.refAlleleCount + s1.refAlleleCount);
            for (const auto& altVal : s0.altAlleleCount)
            {
                depth += altVal.second;
            }
            for (const auto& altVal : s1.altAlleleCount)
            {
                depth += altVal.second;
            }
            iterateMapValue(totalByDepth, depth, obsCount);
        }

        if (key.getStrand0Counts().altAlleleCount.empty() && key.getStrand1Counts().altAlleleCount.empty())
        {
            refOnlyKeyCount++;
        }
        if (key.getStrand0Counts().refAlleleCount==0 && key.getStrand1Counts().refAlleleCount==0)
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
        const SingleSampleContextObservationPattern& key(value.first);
        const auto& obsCount(value.second);
        os << "KEYcount: " << obsCount << "\n";
        os << "KEYS1:\t" << key.getStrand0Counts() << "\n";
        os << "KEYS2:\t" << key.getStrand1Counts() << "\n";
    }
#endif
}



void
ContextData::
merge(
    const ContextData& in)
{
    counts.merge(in.counts);
    excludedRegionSkipped += in.excludedRegionSkipped;
    depthSkipped += in.depthSkipped;
    emptySkipped += in.emptySkipped;
    noiseSkipped += in.noiseSkipped;
}



void
ContextData::
dump(
    std::ostream& os) const
{
    os << "excludedRegionSkippedCount: " << excludedRegionSkipped << "\n";
    os << "depthSkippedCount: " << depthSkipped << "\n";
    os << "emptySkippedCount: " << emptySkipped << "\n";
    os << "noiseSkippedCount: " << noiseSkipped << "\n";

    counts.dump(os);
}



void
Dataset::
addContextInstanceObservation(
    const Context& context,
    const ContextInstanceObservation& observation)
{
    const auto iter(getContextIterator(context));
    iter->second.counts.addSingleSampleContextInstanceObservation(observation);
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
addEmptySkip(
    const Context& context)
{
    const auto iter(getContextIterator(context));
    iter->second.emptySkipped++;
}



void
Dataset::
addNoiseSkip(
    const Context& context)
{
    const auto iter(getContextIterator(context));
    iter->second.noiseSkipped++;
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
    os << "BasecallCounts DUMP_ON\n";
    os << "Total Basecall Contexts: " << _data.size() << "\n";
    for (const auto& value : _data)
    {
        os << "Basecall Context: " << value.first << "\n";
        value.second.dump(os);
    }
    os << "BasecallCounts DUMP_OFF\n";
}



Dataset::data_t::iterator
Dataset::
getContextIterator(
    const Context& context)
{
    const auto insert = _data.insert(std::make_pair(context,data_t::mapped_type()));
    return insert.first;
}

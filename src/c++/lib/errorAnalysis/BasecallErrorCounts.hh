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

#pragma once

#include "boost/serialization/level.hpp"
#include "boost/serialization/map.hpp"

#include <array>
#include <iosfwd>
#include <map>
#include <numeric>
#include <vector>


/// The sequencing context of a basecall, which can be used to segment the basecall error estimation process by context.
///
/// Note that context currently plays a less important role for basecall error analysis compared to that for indels.
struct BasecallErrorContext
{
    BasecallErrorContext() {}

    bool
    operator<(
        const BasecallErrorContext& rhs) const
    {
        return (repeatCount < rhs.repeatCount);
    }

    template<class Archive>
    void serialize(Archive& ar, const unsigned /* version */)
    {
        ar& repeatCount;
    }

    unsigned repeatCount = 1;
};

BOOST_CLASS_IMPLEMENTATION(BasecallErrorContext, boost::serialization::object_serializable)

std::ostream&
operator<<(
    std::ostream& os,
    const BasecallErrorContext& context);


struct StrandBasecallCounts
{
    bool
    operator==(
        const StrandBasecallCounts& rhs) const
    {
        return ((refAlleleCount == rhs.refAlleleCount) && (altAlleleCount == rhs.altAlleleCount));
    }

    bool
    operator<(
        const StrandBasecallCounts& rhs) const
    {
        if (refAlleleCount < rhs.refAlleleCount) return true;
        if (refAlleleCount != rhs.refAlleleCount) return false;
        return (altAlleleCount < rhs.altAlleleCount);
    }

    /// Compress lower bits out of this object's allele counts
    void
    compressCounts();

    typedef std::map<uint16_t,unsigned> qual_count_t;

    /// Reference allele counts
    unsigned refAlleleCount = 0;

    /// Alternate allele counts, stratified by the phred-scaled basecall error probability
    qual_count_t altAlleleCount;
};

std::ostream&
operator<<(std::ostream& os, const StrandBasecallCounts& sbc);


struct BasecallErrorContextObservationData;


/// basecalls are input by strand, but once compressed
/// fwd and rev are rematched to strand0 and strand1 such
/// that strand0 > strand1 -- the per locus strand difference
/// is maintained, but the polarity is thrown away.
///
struct BasecallErrorContextObservation
{
    const StrandBasecallCounts&
    getStrand0Counts() const
    {
        return strand0;
    }

    const StrandBasecallCounts&
    getStrand1Counts() const
    {
        return strand1;
    }

    bool
    operator<(
        const BasecallErrorContextObservation& rhs) const
    {
        if (strand0 < rhs.strand0) return true;
        if (not (strand0 == rhs.strand0)) return false;
        return (strand1 < rhs.strand1);
    }

    void
    compressCounts();

    template<class Archive>
    void serialize(Archive& ar, const unsigned /* version */)
    {
        ar& strand0.refAlleleCount;
        ar& strand0.altAlleleCount;
        ar& strand1.refAlleleCount;
        ar& strand1.altAlleleCount;
    }


private:
    friend BasecallErrorContextObservationData;

    StrandBasecallCounts strand0;
    StrandBasecallCounts strand1;
};

BOOST_CLASS_IMPLEMENTATION(BasecallErrorContextObservation, boost::serialization::object_serializable)


struct BasecallErrorContextObservationData;

/// Basecall count data associated with a single context instance (ie. a single pileup column)
///
/// This is a special version of the data structure used as input into aggregated count structures, where
/// more (sometimes lossy) data compression will be applied.
struct BasecallErrorContextInputObservation
{
    void
    addRefCount(
        const bool isFwdStrand,
        const uint16_t basecallErrorPhredProb);

    void
    addAltCount(
        const bool isFwdStrand,
        const uint16_t basecallErrorPhredProb);

private:
    friend BasecallErrorContextObservationData;

    /// Reference allele observations for forward and reverse strands
    std::array<StrandBasecallCounts::qual_count_t,2> ref;

    /// Alternate allele observations for forward and reverse strands
    std::array<StrandBasecallCounts::qual_count_t,2> alt;
};


/// Basecall error counts associated with a single strand of a single context instance, uncompressed for use by an
/// external error estimation routine
struct BasecallErrorContextObservationExportStrandObservation
{
    bool
    operator==(
        const BasecallErrorContextObservationExportStrandObservation& rhs) const
    {
        return ((refAlleleCount == rhs.refAlleleCount) && (altAlleleCount == rhs.altAlleleCount));
    }

    bool
    operator<(
        const BasecallErrorContextObservationExportStrandObservation& rhs) const
    {
        if (refAlleleCount < rhs.refAlleleCount) return true;
        if (refAlleleCount != rhs.refAlleleCount) return false;
        return (altAlleleCount < rhs.altAlleleCount);
    }

    /// Number of reference allele observations
    unsigned refAlleleCount = 0;

    /// Number of alternate allele observations, stratified by basecall error probability levels (the basecall error
    /// probabilities for each level are stored in
    /// BasecallErrorContextObservationExportData::altAlleleBasecallErrorPhredProbLevels
    std::vector<unsigned> altAlleleCount;
};

/// Basecall error counts associated with a single context instance, uncompressed for use by an
/// external error estimation routine
struct BasecallErrorContextObservationExportObservation
{
    bool
    operator<(
        const BasecallErrorContextObservationExportObservation& rhs) const
    {
        if (strand0 < rhs.strand0) return true;
        if (not (strand0 == rhs.strand0)) return false;
        return (strand1 < rhs.strand1);
    }

    BasecallErrorContextObservationExportStrandObservation strand0;
    BasecallErrorContextObservationExportStrandObservation strand1;
};

/// Basecall error counts for all instances of a single sequence context, uncompressed for use by an
/// external error estimation routine
struct BasecallErrorContextObservationExportData
{
    void
    clear()
    {
        altAlleleBasecallErrorPhredProbLevels.clear();
        refCount.clear();
        observations.clear();
    }

    /// The basecall error levels used for the alternate allele observations
    std::vector<uint16_t> altAlleleBasecallErrorPhredProbLevels;
    /// TODO:
    std::vector<uint64_t> refCount;

    // map value is the number of identical context instances:
    typedef std::map<BasecallErrorContextObservationExportObservation,unsigned> observation_t;
    observation_t observations;
};


struct BasecallErrorData;

/// Basecall counts associated with a specific context
struct BasecallErrorContextObservationData
{
private:
    // map value is the number of identical context instances
    typedef std::map<BasecallErrorContextObservation,unsigned> data_t;
public:
    typedef data_t::const_iterator const_iterator;

    // map key is basecall quality, value is observation counts at that quality level
    typedef std::map<uint16_t, uint64_t> refQual_t;

    void
    addObservation(
        const BasecallErrorContextInputObservation& obs);

    void
    merge(const BasecallErrorContextObservationData& in);

    const_iterator
    begin() const
    {
        return data.begin();
    }

    const_iterator
    end() const
    {
        return data.end();
    }

    const refQual_t&
    getRefQuals() const
    {
        return refAlleleBasecallErrorPhredProbs;
    }

    /// Reformat the data so that it is easier for downstream estimation models to consume it. Primarily this
    /// means undoing various data compression schemes.
    void
    getExportData(BasecallErrorContextObservationExportData& exportData) const;

    void
    dump(std::ostream& os) const;

private:
    friend BasecallErrorData;

    /// Basecall counts stratified by instance count pattern
    data_t data;

    /// Reference allele counts totalled over all instances of the context, but stratified by quality score
    refQual_t refAlleleBasecallErrorPhredProbs;
};

BOOST_CLASS_IMPLEMENTATION(BasecallErrorContextObservationData, boost::serialization::object_serializable)


/// All basecall error data associated with a specific context
struct BasecallErrorData
{
    void
    merge(const BasecallErrorData& in);

    /// debug output
    void
    dump(std::ostream& os) const;

    template<class Archive>
    void serialize(Archive& ar, const unsigned /* version */)
    {
        // adding error.data instead of error here to reduce the total
        // serialization template depth:
        ar& counts.data;
        ar& counts.refAlleleBasecallErrorPhredProbs;
        ar& excludedRegionSkipped;
        ar& depthSkipped;
        ar& emptySkipped;
        ar& noiseSkipped;
    }

    BasecallErrorContextObservationData counts;
    uint64_t excludedRegionSkipped = 0;
    uint64_t depthSkipped = 0;
    uint64_t emptySkipped = 0;
    uint64_t noiseSkipped = 0;
};

BOOST_CLASS_IMPLEMENTATION(BasecallErrorData, boost::serialization::object_serializable)


/// Store all data used for basecall error estimation
struct BasecallErrorCounts
{
private:
    typedef std::map<BasecallErrorContext,BasecallErrorData> data_t;
public:
    typedef data_t::const_iterator const_iterator;

    void
    addSiteObservation(
        const BasecallErrorContext& context,
        const BasecallErrorContextInputObservation& siteObservation);

    void
    addExcludedRegionSkip(
        const BasecallErrorContext& context);

    void
    addDepthSkip(
        const BasecallErrorContext& context);

    void
    addEmptySkip(
        const BasecallErrorContext& context);

    void
    addNoiseSkip(
        const BasecallErrorContext& context);

    void
    merge(const BasecallErrorCounts& in);

    void
    clear()
    {
        _data.clear();
    }

    const_iterator
    begin() const
    {
        return _data.begin();
    }

    const_iterator
    end() const
    {
        return _data.end();
    }

    /// debug output
    void
    dump(std::ostream& os) const;

    template<class Archive>
    void serialize(Archive& ar, const unsigned /* version */)
    {
        ar& _data;
    }

private:
    data_t::iterator
    getContextIterator(
        const BasecallErrorContext& context);

    data_t _data;
};

BOOST_CLASS_IMPLEMENTATION(BasecallErrorCounts, boost::serialization::object_serializable)

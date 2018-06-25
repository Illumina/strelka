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



namespace BasecallCounts
{

/// \brief The sequencing context of a basecall, used to segment the basecall pattern counting and modeling
///
/// Note that context currently plays a less important role for basecall error analysis compared to that for indels.
struct Context
{
    Context() {}

    bool
    operator<(
        const Context& rhs) const
    {
        return (repeatCount < rhs.repeatCount);
    }

    template<class Archive>
    void serialize(
        Archive& ar,
        const unsigned /* version */)
    {
        ar& repeatCount;
    }

    unsigned repeatCount = 1;
};

std::ostream&
operator<<(
    std::ostream& os,
    const Context& context);


/// \brief An allele pattern which could be observed on one strand of any given context instance in one sample
struct SingleSampleSingleStrandContextObservationPattern
{
    bool
    operator==(
        const SingleSampleSingleStrandContextObservationPattern& rhs) const
    {
        return ((refAlleleCount == rhs.refAlleleCount) && (altAlleleCount == rhs.altAlleleCount));
    }

    bool
    operator<(
        const SingleSampleSingleStrandContextObservationPattern& rhs) const
    {
        if (refAlleleCount < rhs.refAlleleCount) return true;
        if (refAlleleCount != rhs.refAlleleCount) return false;
        return (altAlleleCount < rhs.altAlleleCount);
    }

    /// Compress lower bits out of this object's allele counts
    void
    compressCounts();

    typedef std::map<uint16_t, unsigned> qual_count_t;

    /// Reference allele counts
    unsigned refAlleleCount = 0;

    /// Alternate allele counts, stratified by the phred-scaled basecall error probability
    qual_count_t altAlleleCount;
};

std::ostream&
operator<<(
    std::ostream& os,
    const SingleSampleSingleStrandContextObservationPattern& pattern);


struct SingleSampleContextData;

/// \brief An allele pattern which could be observed at any given context instance in one sample
///
/// basecalls are input by strand, but once compressed
/// fwd and rev are rematched to strand0 and strand1 such
/// that strand0 > strand1 -- the per locus strand difference
/// is maintained, but the polarity is thrown away.
///
struct SingleSampleContextObservationPattern
{
    const SingleSampleSingleStrandContextObservationPattern&
    getStrand0Counts() const
    {
        return strand0;
    }

    const SingleSampleSingleStrandContextObservationPattern&
    getStrand1Counts() const
    {
        return strand1;
    }

    bool
    operator<(
        const SingleSampleContextObservationPattern& rhs) const
    {
        if (strand0 < rhs.strand0) return true;
        if (not(strand0 == rhs.strand0)) return false;
        return (strand1 < rhs.strand1);
    }

    void
    compressCounts();

    template<class Archive>
    void serialize(
        Archive& ar,
        const unsigned /* version */)
    {
        ar& strand0.refAlleleCount;
        ar& strand0.altAlleleCount;
        ar& strand1.refAlleleCount;
        ar& strand1.altAlleleCount;
    }


private:
    friend SingleSampleContextData;

    SingleSampleSingleStrandContextObservationPattern strand0;
    SingleSampleSingleStrandContextObservationPattern strand1;
};


struct SingleSampleContextData;

/// \brief Basecall observations associated with a single context instance in the genome
///
/// A single context instance of a basecall context can be thought of as a site in the genome, and the corresponding
/// basecall observations at this site correspond to the data in a single pileup column.
///
/// Basecall counts are simplified to 2 states: reference alleles and non-reference alleles.
///
/// This is a special version of the basecall count data structure used as input into the primary (aggregated) basecall
/// count structures. In the primary structures, a greater level of data compression will be applied.
struct ContextInstanceObservation
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
    friend SingleSampleContextData;

    /// Reference allele observations for forward and reverse strands
    std::array<SingleSampleSingleStrandContextObservationPattern::qual_count_t, 2> ref;

    /// Alternate allele observations for forward and reverse strands
    std::array<SingleSampleSingleStrandContextObservationPattern::qual_count_t, 2> alt;
};


/// \brief An allele pattern which could be observed on one strand of any given context instance in one sample
///
/// The data in this structure are uncompressed for use by an external modeling routine.
struct SingleSampleSingleStrandContextObservationPatternExportFormat
{
    bool
    operator==(
        const SingleSampleSingleStrandContextObservationPatternExportFormat& rhs) const
    {
        return ((refAlleleCount == rhs.refAlleleCount) && (altAlleleCount == rhs.altAlleleCount));
    }

    bool
    operator<(
        const SingleSampleSingleStrandContextObservationPatternExportFormat& rhs) const
    {
        if (refAlleleCount < rhs.refAlleleCount) return true;
        if (refAlleleCount != rhs.refAlleleCount) return false;
        return (altAlleleCount < rhs.altAlleleCount);
    }

    /// Number of reference allele observations
    unsigned refAlleleCount = 0;

    /// Number of alternate allele observations, stratified by basecall error probability levels
    /// (the basecall error probabilities for each level are stored in
    /// BasecallErrorContextObservationExportData::altAlleleBasecallErrorPhredProbLevels)
    std::vector<unsigned> altAlleleCount;
};

/// \brief An allele pattern which could be observed at any given context instance in one sample
///
/// The data in this structure are uncompressed for use by an external modeling routine.
struct SingleSampleContextObservationPatternExportFormat
{
    bool
    operator<(
        const SingleSampleContextObservationPatternExportFormat& rhs) const
    {
        if (strand0 < rhs.strand0) return true;
        if (not(strand0 == rhs.strand0)) return false;
        return (strand1 < rhs.strand1);
    }

    SingleSampleSingleStrandContextObservationPatternExportFormat strand0;
    SingleSampleSingleStrandContextObservationPatternExportFormat strand1;
};

/// \brief All basecall pattern counts associated with a specific context in one sample
///
/// The data in this structure are uncompressed for use by an external modeling routine.
struct SingleSampleContextDataExportFormat
{
    void
    clear()
    {
        altAlleleBasecallErrorPhredProbLevels.clear();
        refCount.clear();
        observations.clear();
    }

    /// The basecall error rate levels used for the alternate allele observations
    std::vector<uint16_t> altAlleleBasecallErrorPhredProbLevels;
    std::vector<uint64_t> refCount;

    // map key is an allele pattern which could be observed at any given context instance in one sample
    // map value is the number of context instances where the pattern described by the map key has been observed
    typedef std::map<SingleSampleContextObservationPatternExportFormat, unsigned> observation_t;
    observation_t observations;
};


struct ContextData;

/// \brief All basecall pattern counts associated with a specific context in one sample
struct SingleSampleContextData
{
private:
    // map key is an allele pattern which could be observed at any given context instance in one sample
    // map value is the number of context instances where the pattern described by the map key has been observed
    typedef std::map<SingleSampleContextObservationPattern, unsigned> data_t;
public:
    typedef data_t::const_iterator const_iterator;

    // map key is basecall quality
    // map value is observation counts at that quality level
    typedef std::map<uint16_t, uint64_t> refQual_t;

    void
    addSingleSampleContextInstanceObservation(
        const ContextInstanceObservation& obs);

    void
    merge(const SingleSampleContextData& in);

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
    exportData(SingleSampleContextDataExportFormat& exportedData) const;

    void
    dump(std::ostream& os) const;

private:
    friend ContextData;

    /// Basecall counts stratified by instance count pattern
    data_t data;

    /// Reference allele counts totalled over all instances of the context, but stratified by quality score
    refQual_t refAlleleBasecallErrorPhredProbs;
};


/// \brief All basecall pattern data associated with a specific context
struct ContextData
{
    void
    merge(const ContextData& in);

    /// debug output
    void
    dump(std::ostream& os) const;

    template<class Archive>
    void serialize(
        Archive& ar,
        const unsigned /* version */)
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

    SingleSampleContextData counts;
    uint64_t excludedRegionSkipped = 0;
    uint64_t depthSkipped = 0;
    uint64_t emptySkipped = 0;
    uint64_t noiseSkipped = 0;
};


/// \brief Store all counts and auxiliary data related to basecalls
struct Dataset
{
private:
    typedef std::map<Context, ContextData> data_t;
public:
    typedef data_t::const_iterator const_iterator;

    void
    addContextInstanceObservation(
        const Context& context,
        const ContextInstanceObservation& observation);

    /// Indicate that a site is skipped because it falls into a user-specified excluded region
    void
    addExcludedRegionSkip(
        const Context& context);

    /// Indicate that a site is skipped due to anomalous depth
    void
    addDepthSkip(
        const Context& context);

    /// Indicate that a site is skipped due to no coverage at the site (after quality filtration)
    void
    addEmptySkip(
        const Context& context);

    /// Indicate that a site is skipped due to excessive sequencing noise in the surrounding locus
    void
    addNoiseSkip(
        const Context& context);

    void
    merge(const Dataset& in);

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
    void serialize(
        Archive& ar,
        const unsigned /* version */)
    {
        ar& _data;
    }

private:
    data_t::iterator
    getContextIterator(
        const Context& context);

    data_t _data;
};

}


// TODO - Prefer to keep these macros inline with each class, work out namespace/macro interactions to allow this.
//
BOOST_CLASS_IMPLEMENTATION(BasecallCounts::Context, boost::serialization::object_serializable)

BOOST_CLASS_IMPLEMENTATION(BasecallCounts::SingleSampleContextObservationPattern, boost::serialization::object_serializable)

BOOST_CLASS_IMPLEMENTATION(BasecallCounts::SingleSampleContextData, boost::serialization::object_serializable)

BOOST_CLASS_IMPLEMENTATION(BasecallCounts::ContextData, boost::serialization::object_serializable)

BOOST_CLASS_IMPLEMENTATION(BasecallCounts::Dataset, boost::serialization::object_serializable)

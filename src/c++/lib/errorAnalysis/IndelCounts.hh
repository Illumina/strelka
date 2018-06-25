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

#include "blt_util/RecordTracker.hh"

#include "boost/serialization/array.hpp"
#include "boost/serialization/level.hpp"
#include "boost/serialization/map.hpp"

#include <array>
#include <iosfwd>
#include <map>
#include <numeric>
#include <vector>
#include <set>


namespace IndelCounts
{

namespace INDEL_TYPE
{
enum index_t
{
    INSERT,
    DELETE,
    SIZE
};

inline
char
symbol(
    const index_t id)
{
    switch (id)
    {
    case INSERT:
        return 'I';
    case DELETE:
        return 'D';
    default:
        return 'X';
    }
}
}


namespace INDEL_SIGNAL_TYPE
{
// second number represents REPEAT_UNITS added or removed
enum index_t
{
    INSERT_1,
    INSERT_2,
    INSERT_GE3,
    DELETE_1,
    DELETE_2,
    DELETE_GE3,
    SIZE
};

inline
const char*
label(
    const index_t id)
{
    switch (id)
    {
    case INSERT_1:
        return "I_1.";
    case INSERT_2:
        return "I_2.";
    case INSERT_GE3:
        return "I_3+";
    case DELETE_1:
        return "D_1.";
    case DELETE_2:
        return "D_2.";
    case DELETE_GE3:
        return "D_3+";
    default:
        return "X_X.";
    }
}

inline
const char*
label(
    const unsigned id)
{
    return label(static_cast<index_t>(id));
}
}


static
GENOTYPE_STATUS::genotype_t
assignStatus(const RecordTracker::indel_value_t& knownVariantOlap)
{
    if (knownVariantOlap.size() == 1)
    {
        // right now, if a known variant has multiple annotations
        // (e.g. OverlapConflict), we will label the variant as
        // UNKNOWN
        return knownVariantOlap.begin()->genotype;
    }

    return GENOTYPE_STATUS::UNKNOWN;
}


/// \brief The sequencing context of an indel, used to segment the indel pattern counting and modeling
struct Context
{
    Context(
        unsigned initRepeatingPatternSize = 1,
        unsigned initRepeatCount = 1)
        : repeatPatternSize(initRepeatingPatternSize),
          repeatCount(initRepeatCount) {}

    bool operator<(const Context& rhs) const
    {
        if (repeatPatternSize < rhs.repeatPatternSize)
        {
            return true;
        }
        if (repeatPatternSize != rhs.repeatPatternSize)
        {
            return false;
        }
        return repeatCount < rhs.repeatCount;
    }

    template<class Archive>
    void serialize(
        Archive& ar,
        const unsigned /* version */)
    {
        ar& repeatPatternSize;
        ar& repeatCount;
    }

    unsigned getRepeatPatternSize() const
    {
        return repeatPatternSize;
    }

    unsigned getRepeatCount() const
    {
        return repeatCount;
    }

private:
    unsigned repeatPatternSize = 1;
    unsigned repeatCount = 1;
};

std::ostream&
operator<<(
    std::ostream& os,
    const Context& context);


/// \brief A non-variant allele pattern which could be observed at any given context instance in one sample
struct SingleSampleNonVariantContextObservationPattern
{
    void
    assignKnownStatus(const RecordTracker::indel_value_t& knownVariantOlap)
    {
        backgroundStatus = assignStatus(knownVariantOlap);
    }

    bool
    operator<(
        const SingleSampleNonVariantContextObservationPattern& rhs) const
    {
        if (depth < rhs.depth) return true;
        if (depth != rhs.depth) return false;
        return (backgroundStatus < rhs.backgroundStatus);
    }

    template<class Archive>
    void serialize(Archive& ar, const unsigned /* version */)
    {
        ar& backgroundStatus& depth;
    }

    unsigned depth;
    GENOTYPE_STATUS::genotype_t backgroundStatus;
};

std::ostream&
operator<<(
    std::ostream& os,
    const SingleSampleNonVariantContextObservationPattern& obs);


struct SingleSampleNonVariantContextData
{
private:
    // map key is a non-variant allele pattern which could be observed at any given context instance in one sample
    // map value is the number of context instances where the pattern described by the map key has been observed
    typedef std::map<SingleSampleNonVariantContextObservationPattern,unsigned> data_t;
public:
    typedef data_t::const_iterator const_iterator;

    void
    addSingleSampleNonVariantInstanceObservation(
        const SingleSampleNonVariantContextObservationPattern& obs);

    void
    merge(const SingleSampleNonVariantContextData& in);

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

    void
    dump(std::ostream& os) const;

    template<class Archive>
    void serialize(Archive& ar, const unsigned /* version */)
    {
        ar& data;
    }

private:
    data_t data;
};


/// \brief Ratio of allele-supporting reads to background depth for each context type
///
/// This ratio is computed for each context and used to estimate the supporting read counts
/// for the reference allele at all non-variant context instances.
struct IndelDepthSupportTotal
{
    IndelDepthSupportTotal(
        unsigned initDepth = 0,
        unsigned initSupport = 0)
        : depth(initDepth), supportCount(initSupport)
    {}

    void
    merge(const IndelDepthSupportTotal& in)
    {
        depth += in.depth;
        supportCount += in.supportCount;
    }

    template<class Archive>
    void serialize(Archive& ar, const unsigned /* version */)
    {
        ar& depth& supportCount;
    }

    double depth;
    double supportCount;
};


/// \brief A candidate variant pattern which could be observed at any given context instance in one sample
struct SingleSampleCandidateVariantContextObservationPattern
{
    SingleSampleCandidateVariantContextObservationPattern()
        : refCount(0)
        , variantStatus(GENOTYPE_STATUS::UNKNOWN)
    {
        std::fill(signalCounts.begin(), signalCounts.end(), 0);
    }

    void
    assignKnownStatus(const RecordTracker::indel_value_t& knownVariantOlap)
    {
        variantStatus = assignStatus(knownVariantOlap);
    }

    unsigned
    totalSignalCount() const
    {
        return std::accumulate(signalCounts.begin(), signalCounts.end(), 0);
    }

    unsigned
    totalCount() const
    {
        return (refCount+totalSignalCount());
    }

    bool
    operator<(
        const SingleSampleCandidateVariantContextObservationPattern& rhs) const
    {
        if (refCount < rhs.refCount) return true;
        if (refCount != rhs.refCount) return false;
        for (unsigned signalIndex(0); signalIndex<INDEL_SIGNAL_TYPE::SIZE; ++signalIndex)
        {
            if (signalCounts[signalIndex] == rhs.signalCounts[signalIndex]) continue;
            return (signalCounts[signalIndex] < rhs.signalCounts[signalIndex]);
        }
        return (variantStatus < rhs.variantStatus);
    }

    template<class Archive>
    void serialize(Archive& ar, const unsigned /* version */)
    {
        ar& variantStatus& refCount& signalCounts;
    }

    /// Number of reads supporting the reference allele
    unsigned refCount;

    /// Optional true genotype label for this locus
    GENOTYPE_STATUS::genotype_t variantStatus;

    /// Number of reads supporting various categories of indel alleles.
    /// The indel allele categories are pre-defined and static.
    std::array<unsigned,INDEL_SIGNAL_TYPE::SIZE> signalCounts;
};

std::ostream&
operator<<(
    std::ostream& os,
    const SingleSampleCandidateVariantContextObservationPattern& obs);


/// \brief All basecall pattern counts associated with a specific context in one sample
///
/// Note that this only includes counts from indel candidates. Counts from non-candidate contexts are recorded
/// in the corresponding "NonVariantContext" data structure.
struct SingleSampleCandidateVariantContextData
{
private:
    // map key is an allele pattern which could be observed at any given context instance in one sample
    // map value is the number of context instances where the pattern described by the map key has been observed
    typedef std::map<SingleSampleCandidateVariantContextObservationPattern, unsigned> data_t;
public:
    typedef data_t::const_iterator const_iterator;

    void
    addSingleSampleCandidateContextInstanceObservation(
        const SingleSampleCandidateVariantContextObservationPattern& obs);

    void
    merge(const SingleSampleCandidateVariantContextData& in);

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

    void
    dump(std::ostream& os) const;

    template<class Archive>
    void serialize(Archive& ar, const unsigned /* version */)
    {
        ar& data;
    }

private:
    data_t data;
};


/// \brief An allele pattern and the count of context instances at which it was observed, for one context in one sample
///
/// The data in this structure are uncompressed for use by an external modeling routine.
struct SingleSampleContextObservationInfoExportFormat
{
    /// The total number of context instances at which the observation pattern described below was observed
    unsigned contextInstanceCount;

    /// Number of supporting observations of the non-alt* allele (*exact definition in flux...)
    double refObservations;

    /// Number of supporting observations of the alt allele
    std::array<unsigned,INDEL_SIGNAL_TYPE::SIZE> altObservations;

    /// Status of indel -- has it been supplied previously as a known variant
    GENOTYPE_STATUS::genotype_t variantStatus = GENOTYPE_STATUS::UNKNOWN;
};

std::ostream&
operator<<(
    std::ostream& os,
    const SingleSampleContextObservationInfoExportFormat& obs);


/// \brief All indel pattern data associated with a specific context
///
/// The data in this structure are uncompressed for use by an external modeling routine.
struct SingleSampleContextDataExportFormat
{
    /// Listing of all observation patterns, and the context instance counts for each pattern.
    ///
    /// Observation patterns are not sorted in any stable order
    std::vector<SingleSampleContextObservationInfoExportFormat> data;
};


/// \brief All indel pattern data associated with a specific context
struct ContextData
{
    /// Processes current data into a format that is more easily
    /// manipulated in a downstream model
    void
    exportData(
        SingleSampleContextDataExportFormat& exportedContextData) const;

    void
    merge(const ContextData& in);

    /// debug output
    void
    dump(std::ostream& os) const;

    template<class Archive>
    void serialize(Archive& ar, const unsigned /* version */)
    {
        ar& nonVariantContextCounts;
        ar& candidateVariantContextCounts;
        ar& depthSupport;
        ar& excludedRegionSkipped;
        ar& depthSkipped;
    }

    SingleSampleNonVariantContextData nonVariantContextCounts;
    SingleSampleCandidateVariantContextData candidateVariantContextCounts;
    IndelDepthSupportTotal depthSupport;
    uint64_t excludedRegionSkipped = 0;
    uint64_t depthSkipped = 0;
};


/// \brief Store all counts and auxiliary data related to indel alleles
struct Dataset
{
private:
    typedef std::map<Context, ContextData> data_t;
public:
    typedef data_t::const_iterator const_iterator;

    /// Submit data for an instance of an indel variant context
    ///
    /// \param depth This is used to estimate read support of the reference allele at non-variant context instances
    void
    addCandidateVariantContextInstanceObservation(
        const Context& context,
        const SingleSampleCandidateVariantContextObservationPattern& candidateVariantContextObservation,
        const unsigned depth);

    /// Submit data for an instance of a non-variant context
    void
    addNonVariantContextInstanceObservation(
        const Context& context,
        const SingleSampleNonVariantContextObservationPattern& nonVariantContextObservation);

    void
    addExcludedRegionSkip(
        const Context& context);

    void
    addDepthSkip(
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

    const_iterator
    find(const Context& context) const
    {
        return _data.find(context);
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
        const Context& context);

    data_t _data;
};

}

// TODO - Prefer to keep these macros inline with each class, work out namespace/macro interactions to allow this.
//
BOOST_CLASS_IMPLEMENTATION(IndelCounts::Context, boost::serialization::object_serializable)

BOOST_CLASS_IMPLEMENTATION(IndelCounts::SingleSampleNonVariantContextObservationPattern, boost::serialization::object_serializable)

BOOST_CLASS_IMPLEMENTATION(IndelCounts::SingleSampleNonVariantContextData, boost::serialization::object_serializable)

BOOST_CLASS_IMPLEMENTATION(IndelCounts::IndelDepthSupportTotal, boost::serialization::object_serializable)

BOOST_CLASS_IMPLEMENTATION(IndelCounts::SingleSampleCandidateVariantContextObservationPattern, boost::serialization::object_serializable)

BOOST_CLASS_IMPLEMENTATION(IndelCounts::SingleSampleCandidateVariantContextData, boost::serialization::object_serializable)

BOOST_CLASS_IMPLEMENTATION(IndelCounts::ContextData, boost::serialization::object_serializable)

BOOST_CLASS_IMPLEMENTATION(IndelCounts::Dataset, boost::serialization::object_serializable)

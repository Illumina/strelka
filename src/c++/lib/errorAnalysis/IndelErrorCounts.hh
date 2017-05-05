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

#pragma once

#include "boost/array.hpp"
#include "boost/serialization/array.hpp"
#include "boost/serialization/level.hpp"
#include "boost/serialization/map.hpp"

#include "blt_util/RecordTracker.hh"

#include <iosfwd>
#include <map>
#include <numeric>
#include <vector>
#include <set>


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

struct IndelErrorContext
{
    IndelErrorContext(
        unsigned initRepeatingPatternSize = 1,
        unsigned initRepeatCount = 1)
        : repeatPatternSize(initRepeatingPatternSize), repeatCount(initRepeatCount)
    {}

    bool operator<(const IndelErrorContext& rhs) const
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
    void serialize(Archive& ar, const unsigned /* version */)
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

BOOST_CLASS_IMPLEMENTATION(IndelErrorContext, boost::serialization::object_serializable)

std::ostream&
operator<<(
    std::ostream& os,
    const IndelErrorContext& context);


struct IndelBackgroundObservation
{
    void
    assignKnownStatus(const RecordTracker::indel_value_t& knownVariantOlap)
    {
        backgroundStatus = assignStatus(knownVariantOlap);
    }

    bool
    operator<(
        const IndelBackgroundObservation& rhs) const
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
    const IndelBackgroundObservation& obs);

BOOST_CLASS_IMPLEMENTATION(IndelBackgroundObservation, boost::serialization::object_serializable)


struct IndelBackgroundObservationData
{
private:
    typedef std::map<IndelBackgroundObservation,unsigned> data_t;
public:
    typedef data_t::const_iterator const_iterator;

    void
    addObservation(
        const IndelBackgroundObservation& obs);

    void
    merge(const IndelBackgroundObservationData& in);

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
    // value is number of observations:
    std::map<IndelBackgroundObservation,unsigned> data;
};

BOOST_CLASS_IMPLEMENTATION(IndelBackgroundObservationData, boost::serialization::object_serializable)


/// this object is used to construct a single support/depth ratio for each context type
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

BOOST_CLASS_IMPLEMENTATION(IndelDepthSupportTotal, boost::serialization::object_serializable)


struct IndelErrorContextObservation
{
    IndelErrorContextObservation()
        : refCount(0)
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
        const IndelErrorContextObservation& rhs) const
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

    unsigned refCount;
    GENOTYPE_STATUS::genotype_t variantStatus;
    /// note: can't get std::array to serialize correctly on clang, so using boost::array instead
    boost::array<unsigned,INDEL_SIGNAL_TYPE::SIZE> signalCounts;
};

std::ostream&
operator<<(
    std::ostream& os,
    const IndelErrorContextObservation& obs);

BOOST_CLASS_IMPLEMENTATION(IndelErrorContextObservation, boost::serialization::object_serializable)


struct IndelErrorContextObservationData
{
private:
    typedef std::map<IndelErrorContextObservation,unsigned> data_t;
public:
    typedef data_t::const_iterator const_iterator;

    void
    addObservation(
        const IndelErrorContextObservation& obs);

    void
    merge(const IndelErrorContextObservationData& in);

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
    // value is number of observations:
    std::map<IndelErrorContextObservation,unsigned> data;
};

BOOST_CLASS_IMPLEMENTATION(IndelErrorContextObservationData, boost::serialization::object_serializable)



/// this struct is used for output to downstream models only
///
struct ExportedIndelObservations
{
    /// the number of times we observe the observation pattern
    unsigned observationCount;

    /// number of supporting observations of the non-alt* allele (*exact definition in flux...)
    double refObservations;

    /// number of supporting observations of the alt allele
    boost::array<unsigned,INDEL_SIGNAL_TYPE::SIZE> altObservations;

    /// status of indel -- has it been supplied previously as a known variant
    GENOTYPE_STATUS::genotype_t variantStatus = GENOTYPE_STATUS::UNKNOWN;
};

std::ostream&
operator<<(
    std::ostream& os,
    const ExportedIndelObservations& obs);


struct IndelErrorData
{
    /// processes current data into a format that is more easily
    /// manipulated in a downstream model
    void
    exportObservations(
        std::vector<ExportedIndelObservations>& counts) const;
    void
    merge(const IndelErrorData& in);

    /// debug output
    void
    dump(std::ostream& os) const;

    template<class Archive>
    void serialize(Archive& ar, const unsigned /* version */)
    {
        ar& background;
        ar& error;
        ar& depthSupport;
        ar& excludedRegionSkipped;
        ar& depthSkipped;
    }

    IndelBackgroundObservationData background;
    IndelErrorContextObservationData error;
    IndelDepthSupportTotal depthSupport;
    uint64_t excludedRegionSkipped = 0;
    uint64_t depthSkipped = 0;
};

BOOST_CLASS_IMPLEMENTATION(IndelErrorData, boost::serialization::object_serializable)



struct IndelErrorCounts
{
private:
    typedef std::map<IndelErrorContext,IndelErrorData> data_t;
public:
    typedef data_t::const_iterator const_iterator;

    void
    addError(
        const IndelErrorContext& context,
        const IndelErrorContextObservation& errorObservation,
        const unsigned depth);

    void
    addBackground(
        const IndelErrorContext& context,
        const IndelBackgroundObservation& backgroundObservation);

    void
    addExcludedRegionSkip(
        const IndelErrorContext& context);

    void
    addDepthSkip(
        const IndelErrorContext& context);

    void
    merge(const IndelErrorCounts& in);

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
    find(const IndelErrorContext& key) const
    {
        return _data.find(key);
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
        const IndelErrorContext& context);

    data_t _data;
};

BOOST_CLASS_IMPLEMENTATION(IndelErrorCounts, boost::serialization::object_serializable)

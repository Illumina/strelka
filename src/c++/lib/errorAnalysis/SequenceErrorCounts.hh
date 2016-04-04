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

#pragma once

#include "boost/serialization/level.hpp"
#include "boost/array.hpp"

#include <iosfwd>
#include <map>
#include <numeric>
#include <vector>


namespace INDEL_TYPE
{
    enum index_t {
        INSERT,
        DELETE,
        SIZE
    };

    inline
    char
    symbol(
        const index_t id)
    {
        switch(id)
        {
        case INSERT: return 'I';
        case DELETE: return 'D';
        default:     return 'X';
        }
    }
}


namespace SIGNAL_TYPE
{
    // second number represents REPEAT_UNITS added or removed
    enum index_t {
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
        switch(id)
        {
        case INSERT_1:   return "I_1.";
        case INSERT_2:   return "I_2.";
        case INSERT_GE3: return "I_3+";
        case DELETE_1:   return "D_1.";
        case DELETE_2:   return "D_2.";
        case DELETE_GE3: return "D_3+";
        default:     return "X_X.";
        }
    }
}


struct SequenceErrorContext
{
    bool
    operator<(
        const SequenceErrorContext& rhs) const
    {
        return (repeatCount < rhs.repeatCount);
    }

    template<class Archive>
    void serialize(Archive& ar, const unsigned /* version */)
    {
        ar& repeatCount;
    }

    unsigned repeatCount;
};

BOOST_CLASS_IMPLEMENTATION(SequenceErrorContext, boost::serialization::object_serializable)

std::ostream&
operator<<(
    std::ostream& os,
    const SequenceErrorContext& context);


struct SequenceBackgroundObservation
{
    bool
    operator<(
        const SequenceBackgroundObservation& rhs) const
    {
        return (depth < rhs.depth);
    }

    template<class Archive>
    void serialize(Archive& ar, const unsigned /* version */)
    {
        ar& depth;
    }

    unsigned depth;
};

BOOST_CLASS_IMPLEMENTATION(SequenceBackgroundObservation, boost::serialization::object_serializable)


struct SequenceBackgroundObservationData
{
private:
    typedef std::map<SequenceBackgroundObservation,unsigned> data_t;
public:
    typedef data_t::const_iterator const_iterator;

    void
    addObservation(
        const SequenceBackgroundObservation& obs);

    void
    merge(const SequenceBackgroundObservationData& in);

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
    std::map<SequenceBackgroundObservation,unsigned> data;
};

BOOST_CLASS_IMPLEMENTATION(SequenceBackgroundObservationData, boost::serialization::object_serializable)


/// this object is used to construct a single support/depth ratio for each context type
struct SequenceDepthSupportTotal
{
    SequenceDepthSupportTotal(
        unsigned initDepth = 0,
        unsigned initSupport = 0)
        : depth(initDepth), supportCount(initSupport)
    {}

    void
    merge(const SequenceDepthSupportTotal& in)
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

BOOST_CLASS_IMPLEMENTATION(SequenceDepthSupportTotal, boost::serialization::object_serializable)


struct SequenceErrorContextObservation
{
    SequenceErrorContextObservation()
      : refCount(0)
    {
        std::fill(signalCounts.begin(), signalCounts.end(), 0);
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
        const SequenceErrorContextObservation& rhs) const
    {
        if (refCount == rhs.refCount)
        {
            for (unsigned signalIndex(0); signalIndex<SIGNAL_TYPE::SIZE; ++signalIndex)
            {
                if (signalCounts[signalIndex] == rhs.signalCounts[signalIndex]) continue;
                return (signalCounts[signalIndex] < rhs.signalCounts[signalIndex]);
            }
        }
        return (refCount < rhs.refCount);
    }

    template<class Archive>
    void serialize(Archive& ar, const unsigned /* version */)
    {
        ar& refCount& signalCounts;
    }

    unsigned refCount;
    /// note: can't get std::array to serialize correctly on clang
    boost::array<unsigned,SIGNAL_TYPE::SIZE> signalCounts;
};

BOOST_CLASS_IMPLEMENTATION(SequenceErrorContextObservation, boost::serialization::object_serializable)


struct SequenceErrorContextObservationData
{
private:
    typedef std::map<SequenceErrorContextObservation,unsigned> data_t;
public:
    typedef data_t::const_iterator const_iterator;

    void
    addObservation(
        const SequenceErrorContextObservation& obs);

    void
    merge(const SequenceErrorContextObservationData& in);

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
    std::map<SequenceErrorContextObservation,unsigned> data;
};

BOOST_CLASS_IMPLEMENTATION(SequenceErrorContextObservationData, boost::serialization::object_serializable)



/// this struct is used for output to downstream models only
///
struct ExportedObservations
{
    /// the number of times we observe the observation pattern below:
    unsigned repeatCount;

    /// number of supporting observations of the non-alt* allele (*exact definition in flux...)
    unsigned refObservations;

    /// number of supporting observations of the alt allele
    boost::array<unsigned,SIGNAL_TYPE::SIZE> altObservations;
};


struct SequenceErrorData
{
    /// processes current data into a format that is more easily
    /// manipulated in a downstream model
    void
    exportObservations(
        std::vector<ExportedObservations>& counts) const;

    void
    merge(const SequenceErrorData& in);

    /// debug output
    void
    dump(std::ostream& os) const;

    template<class Archive>
    void serialize(Archive& ar, const unsigned /* version */)
    {
        ar& background& error& depthSupport& skipped;
    }

    SequenceBackgroundObservationData background;
    SequenceErrorContextObservationData error;
    SequenceDepthSupportTotal depthSupport;
    unsigned skipped = 0;
};

BOOST_CLASS_IMPLEMENTATION(SequenceErrorData, boost::serialization::object_serializable)



/// Stores sequencing error patterns observed from data
///
/// Used by downstream operations to estimate useful error
/// parameters for modeling, but no parameters are directly
/// stored in this object, this object is only a compressed
/// abstraction of the input data
///
struct SequenceErrorCounts
{
private:
    typedef std::map<SequenceErrorContext,SequenceErrorData> data_t;
public:
    typedef data_t::const_iterator const_iterator;

    void
    addError(
        const SequenceErrorContext& context,
        const SequenceErrorContextObservation& errorObservation,
        const unsigned depth);

    void
    addBackground(
        const SequenceErrorContext& context,
        const SequenceBackgroundObservation& backgroundObservation);

    void
    addDepthSkip(
        const SequenceErrorContext& context);

    void
    merge(const SequenceErrorCounts& in);

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

    void
    save(const char* filename) const;

    void
    load(const char* filename);

    /// debug output
    void
    dump(std::ostream& os) const;

private:
    data_t::iterator
    getContextIterator(
        const SequenceErrorContext& context);

    data_t _data;
};

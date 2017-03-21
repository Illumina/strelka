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

#include "boost/serialization/level.hpp"
#include "boost/serialization/map.hpp"

#include <array>
#include <iosfwd>
#include <map>
#include <numeric>
#include <vector>


struct BaseErrorContext
{
    BaseErrorContext() {}

    bool
    operator<(
        const BaseErrorContext& rhs) const
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

BOOST_CLASS_IMPLEMENTATION(BaseErrorContext, boost::serialization::object_serializable)

std::ostream&
operator<<(
    std::ostream& os,
    const BaseErrorContext& context);


struct StrandBaseCounts
{
    bool
    operator==(
        const StrandBaseCounts& rhs) const
    {
        return ((refCount == rhs.refCount) && (alt == rhs.alt));
    }

    bool
    operator<(
        const StrandBaseCounts& rhs) const
    {
        if (refCount < rhs.refCount) return true;
        if (refCount != rhs.refCount) return false;
        return (alt < rhs.alt);
    }

    void
    compressCounts();

    typedef std::map<uint16_t,unsigned> qual_count_t;

    unsigned refCount = 0;
    qual_count_t alt;
};

std::ostream&
operator<<(std::ostream& os, const StrandBaseCounts& sbc);


struct BaseErrorContextObservationData;


/// basecalls are input by strand, but once compressed
/// fwd and rev are rematched to strand0 and strand1 such
/// that strand0 > strand1 -- the per locus strand difference
/// is maintained, but the polarity is thrown away.
///
struct BaseErrorContextObservation
{
    const StrandBaseCounts&
    getStrand0Counts() const
    {
        return strand0;
    }

    const StrandBaseCounts&
    getStrand1Counts() const
    {
        return strand1;
    }

    bool
    operator<(
        const BaseErrorContextObservation& rhs) const
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
        ar& strand0.refCount;
        ar& strand0.alt;
        ar& strand1.refCount;
        ar& strand1.alt;
    }


private:
    friend BaseErrorContextObservationData;

    StrandBaseCounts strand0;
    StrandBaseCounts strand1;
};

BOOST_CLASS_IMPLEMENTATION(BaseErrorContextObservation, boost::serialization::object_serializable)


struct BaseErrorContextObservationData;

/// special version of struct only used for input from client code:
struct BaseErrorContextInputObservation
{
    void
    addRefCount(
        const bool isFwdStrand,
        const uint16_t qual);

    void
    addAltCount(
        const bool isFwdStrand,
        const uint16_t qual);

private:
    friend BaseErrorContextObservationData;
    std::array<StrandBaseCounts::qual_count_t,2> ref;
    std::array<StrandBaseCounts::qual_count_t,2> alt;
};


///
/// export structs, designed to provide a simpler analysis alternative for client
/// code only:
///
struct BaseErrorContextObservationExportStrandObservation
{
    bool
    operator==(
        const BaseErrorContextObservationExportStrandObservation& rhs) const
    {
        return ((refCount == rhs.refCount) && (altCount == rhs.altCount));
    }

    bool
    operator<(
        const BaseErrorContextObservationExportStrandObservation& rhs) const
    {
        if (refCount < rhs.refCount) return true;
        if (refCount != rhs.refCount) return false;
        return (altCount < rhs.altCount);
    }

    unsigned refCount = 0;
    std::vector<unsigned> altCount;
};

struct BaseErrorContextObservationExportObservation
{
    bool
    operator<(
        const BaseErrorContextObservationExportObservation& rhs) const
    {
        if (strand0 < rhs.strand0) return true;
        if (not (strand0 == rhs.strand0)) return false;
        return (strand1 < rhs.strand1);
    }

    BaseErrorContextObservationExportStrandObservation strand0;
    BaseErrorContextObservationExportStrandObservation strand1;
};

struct BaseErrorContextObservationExportData
{
    void
    clear()
    {
        qualLevels.clear();
        refCount.clear();
        observations.clear();
    }

    std::vector<uint16_t> qualLevels;
    std::vector<uint64_t> refCount;

    typedef std::map<BaseErrorContextObservationExportObservation,unsigned> observation_t;
    observation_t observations;

};


struct BaseErrorData;

struct BaseErrorContextObservationData
{
private:
    typedef std::map<BaseErrorContextObservation,unsigned> data_t;
public:
    typedef data_t::const_iterator const_iterator;

    typedef std::map<uint16_t, uint64_t> refQual_t;

    void
    addObservation(
        const BaseErrorContextInputObservation& obs);

    void
    merge(const BaseErrorContextObservationData& in);

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
        return refQuals;
    }

    // reformat the results in a way that's easier for downstream models to consume:
    void
    getExportData(BaseErrorContextObservationExportData& exportData) const;

    void
    dump(std::ostream& os) const;

private:
    friend BaseErrorData;

    // value is number of observations:
    data_t data;
    refQual_t refQuals;
};

BOOST_CLASS_IMPLEMENTATION(BaseErrorContextObservationData, boost::serialization::object_serializable)


struct BaseErrorData
{
    void
    merge(const BaseErrorData& in);

    /// debug output
    void
    dump(std::ostream& os) const;

    template<class Archive>
    void serialize(Archive& ar, const unsigned /* version */)
    {
        // adding error.data and not error here to reduce total
        // serialize template depth:
        ar& error.data;
        ar& error.refQuals;
        ar& excludedRegionSkipped;
        ar& depthSkipped;
        ar& emptySkipped;
        ar& noiseSkipped;
    }

    BaseErrorContextObservationData error;
    uint64_t excludedRegionSkipped = 0;
    uint64_t depthSkipped = 0;
    uint64_t emptySkipped = 0;
    uint64_t noiseSkipped = 0;
};

BOOST_CLASS_IMPLEMENTATION(BaseErrorData, boost::serialization::object_serializable)



struct BaseErrorCounts
{
private:
    typedef std::map<BaseErrorContext,BaseErrorData> data_t;
public:
    typedef data_t::const_iterator const_iterator;

    void
    addSiteObservation(
        const BaseErrorContext& context,
        const BaseErrorContextInputObservation& siteObservation);

    void
    addExcludedRegionSkip(
        const BaseErrorContext& context);

    void
    addDepthSkip(
        const BaseErrorContext& context);

    void
    addEmptySkip(
        const BaseErrorContext& context);

    void
    addNoiseSkip(
        const BaseErrorContext& context);

    void
    merge(const BaseErrorCounts& in);

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
        const BaseErrorContext& context);

    data_t _data;
};

BOOST_CLASS_IMPLEMENTATION(BaseErrorCounts, boost::serialization::object_serializable)

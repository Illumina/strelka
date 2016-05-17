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
/// \author Mitch Bekritsky
///

#pragma once

#include "blt_util/known_pos_range2.hh"
#include "blt_util/RegionTracker.hh"
#include "blt_util/blt_types.hh"

#include "htsapi/vcf_record.hh"

#include "boost/icl/discrete_interval.hpp"
#include "boost/icl/interval_map.hpp"

#include <iosfwd>
#include <set>

struct RecordTracker
{
    bool
    empty() const
    {
        return _records.empty();
    }

    void
    clear()
    {
        _records.clear();
    }

    bool
    intersectingRecord(
        const pos_t pos,
        std::set<std::string>& records) const
    {
        return intersectingRecordImpl(pos, pos + 1, records);
    }

    bool
    intersectingRecord(
        const known_pos_range2 range,
        std::set<std::string>& records) const
    {
        return intersectingRecordImpl(range.begin_pos(), range.end_pos(), records);
    }

    void
    addVcfRecord(
        const vcf_record& vcfRecord);

    void 
    dump(
        std::ostream& os) const;

    unsigned
    size() const
    {
        return _records.size();
    }

    typedef std::set<std::string> value_t;

private:
    typedef boost::icl::interval_map<pos_t, std::set<std::string> > record_t;
    typedef boost::icl::discrete_interval<pos_t> interval_t;

    bool
    intersectingRecordImpl(
        const pos_t beginPos,
        const pos_t endPos,
        std::set<std::string>& records) const;

    record_t _records;
};

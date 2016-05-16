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

#include <set>
#include <iosfwd>

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
        std::string& record) const
    {
        return intersectingRecordImpl(pos, pos + 1, record);
    }

    bool
    intersectingRecord(
        const known_pos_range2 range,
        std::string& record) const
    {
        return intersectingRecordImpl(range.begin_pos(), range.end_pos(), record);
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

    typedef std::map<known_pos_range2, std::string, PosRangeEndSort> record_t;

private:

    bool
    intersectingRecordImpl(
        const pos_t beginPos,
        const pos_t endPos,
        std::string& record) const;

    record_t _records;
};

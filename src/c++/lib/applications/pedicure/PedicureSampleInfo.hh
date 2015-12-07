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

#include <cassert>

#include <vector>


namespace PEDICURE_SAMPLETYPE
{
enum index_t
{
    PROBAND,
    PARENT,
    SIBLING,
    SIZE
};

inline
const char*
get_label(const unsigned i)
{
    switch (static_cast<index_t>(i))
    {
    case PROBAND:
        return "PROBAND";
    case PARENT:
        return "PARENT";
    case SIBLING:
        return "SIBLING";
    default:
        return "UNKNOWN";
    }
}
}

/// possible future expansion of type def:
namespace PEDICURE_SEX
{
enum index_t
{
    UNKNOWN,
    MALE,
    FEMALE,
    SIZE
};
}


/// tracks all sample information provided by the user
struct SampleInfo
{
    /// relationship of sample to proband there's no use for an unknown value here:
    PEDICURE_SAMPLETYPE::index_t stype = PEDICURE_SAMPLETYPE::PROBAND;
};



struct SampleInfoManager
{
private:
    typedef std::vector<SampleInfo> dtype;
public:
    typedef dtype::const_iterator const_iterator;

    unsigned
    size() const
    {
        return _data.size();
    }

    const_iterator
    begin() const
    {
        return _data.cbegin();
    }

    const_iterator
    end() const
    {
        return _data.cend();
    }

    const std::vector<unsigned>&
    getTypeIndexList(const PEDICURE_SAMPLETYPE::index_t i) const
    {
        assert((i >= 0) && (i < PEDICURE_SAMPLETYPE::SIZE));
        return _sindex[i];
    }

    const SampleInfo&
    getSampleInfo(const unsigned index) const
    {
        assert(index < size());
        return _data[index];
    }

    void
    push_back(const SampleInfo& info)
    {
        _sindex[info.stype].push_back(_data.size());
        _data.push_back(info);
    }

private:
    typedef std::vector<unsigned> stypeIndex;
    stypeIndex _sindex[PEDICURE_SAMPLETYPE::SIZE];

    dtype _data;
};

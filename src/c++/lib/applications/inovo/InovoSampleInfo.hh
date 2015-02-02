// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Starka
// Copyright (c) 2009-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

///
/// \author Chris Saunders
///

#pragma once

#include <cassert>


namespace INOVO_SAMPLETYPE
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

namespace INOVO_GENDER
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
    INOVO_SAMPLETYPE::index_t stype = INOVO_SAMPLETYPE::PROBAND;

    /// gender of sample, this is provided for future expansions but will be ignored for POC
    INOVO_GENDER::index_t gtype = INOVO_GENDER::UNKNOWN;
};



struct SampleInfoManager
{
    unsigned
    size() const
    {
        return _data.size();
    }

    const std::vector<unsigned>&
    getTypeIndexList(const INOVO_SAMPLETYPE::index_t i)
    {
        assert((i >= 0) && (i < INOVO_SAMPLETYPE::SIZE));
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
    stypeIndex _sindex[INOVO_SAMPLETYPE::SIZE];

    std::vector<SampleInfo> _data;
};

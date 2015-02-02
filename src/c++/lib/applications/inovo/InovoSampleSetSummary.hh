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

/// \file
///
/// \author Chris Saunders
///

///
/// note coding convention for all ranges '_pos fields' is:
/// XXX_begin_pos is zero-indexed position at the beginning of the range
/// XXX_end_pos is zero-index position 1 step after the end of the range
///

#pragma once

#include "inovo_shared.hh"
#include "starling_common/SampleSetSummary.hh"


struct InovoSampleSetSummary : public SampleSetSummary
{
    InovoSampleSetSummary(
        const inovo_options& opt)
    : SampleSetSummary(),
      _sinfo(opt.alignFileOpt.alignmentSampleInfo)
    {}

    unsigned
    size() const override
    {
        return _sinfo.size();
    }

    const char*
    get_label(
        const unsigned i) const override
    {
        assert(i<size());
        return INOVO_SAMPLETYPE::get_label(_sinfo.getSampleInfo(i).stype);
    }

    const char*
    get_prefix(
        const unsigned i,
        const bool /*is_tier1*/) const override
    {
        assert(i<size());
        return get_label(i);
    }
private:
    const SampleInfoManager& _sinfo;
};

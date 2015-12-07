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

#include "pedicure_shared.hh"
#include "starling_common/SampleSetSummary.hh"


struct PedicureSampleSetSummary : public SampleSetSummary
{
    explicit
    PedicureSampleSetSummary(
        const pedicure_options& opt)
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
        return PEDICURE_SAMPLETYPE::get_label(_sinfo.getSampleInfo(i).stype);
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

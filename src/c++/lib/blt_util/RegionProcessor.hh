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

#include "blt_util/blt_types.hh"
#include "blt_util/pos_range.hh"

#include <iosfwd>
#include <string>


/// manage creation of a region track written out in BED format
struct RegionProcessor
{
    explicit
    RegionProcessor(
        std::ostream* osptr) :
        _osptr(osptr)
    {}

    ~RegionProcessor()
    {
        flush();
    }

    /// positions must be added in order:
    void
    addToRegion(
        const std::string& chrom,
        const pos_t outputPos);

    // write out any pending ranges:
    void
    flush();

private:
    std::ostream* _osptr;

    bool _is_range = false;
    std::string _chrom;
    pos_range _prange;
};

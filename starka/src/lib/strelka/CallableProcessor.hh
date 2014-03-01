// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2013 Illumina, Inc.
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

#pragma once

#include "blt_util/blt_types.hh"
#include "blt_util/pos_range.hh"
#include "strelka/position_somatic_snv_grid_shared.hh"

#include <iosfwd>
#include <string>


/// manage creation of a callable bed track
struct CallableProcessor
{
    CallableProcessor(
        std::ostream* osptr = NULL) :
        _minQSS(15),
        _minNQSS(25),
        _osptr(osptr),
        _is_range(false)
    {}

    ~CallableProcessor()
    {
        flush();
    }

    void
    add(
        const std::string& chrom,
        const pos_t outputPos,
        const somatic_snv_genotype_grid& sgtg);

    // write out any pending ranges:
    void
    flush();

private:
    const int _minQSS;
    const int _minNQSS;
    std::ostream* _osptr;

    bool _is_range;
    std::string _chrom;
    pos_range _prange;
};

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

#include "position_somatic_snv_grid_shared.hh"

#include "blt_util/RegionProcessor.hh"


/// manage creation of a somatic callable bed track
struct SomaticCallableProcessor : public RegionProcessor
{
    typedef RegionProcessor base_t;

    explicit
    SomaticCallableProcessor(
        std::ostream* osptr) :
        base_t(osptr)
    {}

    void
    addToRegion(
        const std::string& chrom,
        const pos_t outputPos,
        const somatic_snv_genotype_grid& sgtg)
    {
        if ((sgtg.rs.snv_qphred < _minQSS) && (sgtg.rs.nonsomatic_qphred < _minNQSS)) return;
        base_t::addToRegion(chrom,outputPos);
    }

private:
    static const int _minQSS = 15;
    static const int _minNQSS = 15;
};

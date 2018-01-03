//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2018 Illumina, Inc.
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

#include "somatic_result_set.hh"

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
        if ((sgtg.rs.qphred < _minQSS) && (sgtg.rs.nonsomatic_qphred < _minNQSS)) return;
        base_t::addToRegion(chrom,outputPos);
    }

private:
    static const int _minQSS = 15;
    static const int _minNQSS = 15;
};

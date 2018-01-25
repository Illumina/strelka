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

#pragma once

#include "starling_shared.hh"
#include "blt_util/chrom_depth_map.hh"

#include <iosfwd>


/// Add additional content to the VCF header specific to germline analysis outputs
///
/// \param isGenomeVCF If true, the header output is formatted for gVCF. Otherwise the header output is formatted
///                    for a conventional variants VCF file.
///
void
finishGermlineVCFheader(
    const starling_options& opt,
    const gvcf_deriv_options& dopt,
    const cdmap_t& chrom_depth,
    const std::vector<std::string>& sampleNames,
    const bool isGenomeVCF,
    std::ostream& os);

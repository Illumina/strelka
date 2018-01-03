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

#include "vcf_streamer.hh"
#include "blt_util/reference_contig_segment.hh"


/// Return true if the vcf record REF field is compatible with the given reference
///
/// "Compatible with the reference" means that each base in the VCF REF field either matches that in the reference, or
/// the VCF REF base is 'N' or the reference base is 'N'.
bool
isExpectedVcfReference(
    const reference_contig_segment& ref,
    const vcf_record& vcfRecord);

/// Test if the vcf record REF field matches the given reference and throw if it does not
///
/// \param[in] vcfStreamer Provides the current vcf_record and additional context in exception message
void
assertExpectedVcfReference(
    const reference_contig_segment& ref,
    const vcf_streamer& vcfStreamer);

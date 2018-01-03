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

/// \file
/// \author Chris Saunders
///

#pragma once

#include "htsapi/bam_seq.hh"

/// Return the number of contiguous N's from the late-cycle edge of the read
///
/// For example, for read="NNACGTNNN",
/// ...the return value is 3 if isFwdStrand=true and 2 if ifFwdStrand=false
unsigned
getReadAmbiguousEndLength(
    const bam_seq& bseq,
    const bool isFwdStrand);

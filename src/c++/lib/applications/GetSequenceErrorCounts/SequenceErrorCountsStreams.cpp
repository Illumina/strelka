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

/// \author Chris Saunders
///

#include "SequenceErrorCountsStreams.hh"
#include "htsapi/bam_header_util.hh"

#include <cassert>

#include <fstream>
#include <iostream>



SequenceErrorCountsStreams::
SequenceErrorCountsStreams(
    const SequenceErrorCountsOptions& opt,
    const prog_info& pinfo,
    const bam_hdr_t& header,
    const SampleSetSummary& ssi)
    : base_t(opt,pinfo,ssi)
{
    assert(_n_samples == 1);
    _sampleName = get_bam_header_sample_name(header);
}

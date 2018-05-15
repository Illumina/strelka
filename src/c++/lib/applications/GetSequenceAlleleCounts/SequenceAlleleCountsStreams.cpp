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

#include "SequenceAlleleCountsStreams.hh"
#include "htsapi/bam_header_util.hh"

#include <cassert>

#include <fstream>
#include <iostream>



SequenceAlleleCountsStreams::
SequenceAlleleCountsStreams(
    const SequenceAlleleCountsOptions& opt,
    const prog_info& pinfo,
    const bam_hdr_t& header)
    : base_t(1)
{
    assert(getSampleCount() == 1);
    _sampleName = get_bam_header_sample_name(header);

    if (opt.is_write_observations())
    {
        std::ofstream* fosptr(new std::ofstream);
        _observation_bed_osptr.reset(fosptr);
        std::ofstream& fos(*fosptr);
        open_ofstream(pinfo,opt.observationsBedFilename,"obs_bed",fos);
    }
}

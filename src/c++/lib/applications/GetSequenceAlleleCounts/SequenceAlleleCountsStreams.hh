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

#include "SequenceAlleleCountsOptions.hh"
#include "starling_common/starling_streams_base.hh"


struct SequenceAlleleCountsStreams : public starling_streams_base
{
    typedef starling_streams_base base_t;

    SequenceAlleleCountsStreams(
        const SequenceAlleleCountsOptions& client_opt,
        const prog_info& pinfo,
        const bam_hdr_t& bam_header);

    const std::string&
    getSampleName() const
    {
        return _sampleName;
    }

    std::ostream*
    observation_bed_osptr() const
    {
        return _observation_bed_osptr.get();
    }

private:
    std::string _sampleName;
    std::unique_ptr<std::ostream> _observation_bed_osptr;
};

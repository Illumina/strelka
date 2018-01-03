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

/// \author Chris Saunders
///

#pragma once

#include "snoise_shared.hh"
#include "starling_common/starling_streams_base.hh"


struct snoise_streams : public starling_streams_base
{
    typedef starling_streams_base base_t;

    snoise_streams(
        const snoise_options& client_opt,
        const prog_info& pinfo,
        const bam_hdr_t& bam_header,
        const unsigned sampleCount);

    std::ostream*
    snoise_osptr() const
    {
        return _snoise_osptr;
    }

    std::ostream* _snoise_osptr;
};

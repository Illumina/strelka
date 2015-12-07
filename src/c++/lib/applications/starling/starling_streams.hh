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

#pragma once

#include "starling_shared.hh"
#include "starling_common/starling_streams_base.hh"


struct starling_streams : public starling_streams_base
{
    typedef starling_streams_base base_t;

    starling_streams(const starling_options& client_opt,
                     const prog_info& pinfo,
                     const bam_header_t* const bam_header,
                     const SampleSetSummary& ssi);

    std::ostream*
    gvcf_osptr() const
    {
        return _gvcf_osptr;
    }

private:
    static
    std::ostream*
    initialize_gvcf_file(
        const starling_options& opt,
        const prog_info& pinfo,
        const std::string& filename,
        const bam_header_t* const header,
        std::unique_ptr<std::ostream>& os_ptr_auto);

    std::ostream* _gvcf_osptr;
    std::unique_ptr<std::ostream> _gvcf_osptr_auto;
};

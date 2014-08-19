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

/// \author Chris Saunders
///

#pragma once

#include "starling_common/starling_streams_base.hh"


struct snoise_streams : public starling_streams_base
{
    typedef starling_streams_base base_t;

    snoise_streams(
        const snoise_options& client_opt,
        const prog_info& pinfo,
        const bam_header_t* const bam_header);

    std::ostream*
    snoise_osptr() const
    {
        return _snoise_osptr.get();
    }

    std::unique_ptr<std::ostream> _snoise_osptr;
};

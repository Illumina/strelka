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

#include "inovo_shared.hh"
#include "InovoSampleSetSummary.hh"

#include "starling_common/starling_streams_base.hh"



struct inovo_streams : public starling_streams_base
{
    typedef starling_streams_base base_t;

    inovo_streams(
        const inovo_options& opt,
        const inovo_deriv_options& dopt,
        const prog_info& pinfo,
        const bam_header_t* const bam_header,
        const InovoSampleSetSummary& ssi);

    std::ostream*
    denovo_osptr() const
    {
        return _denovo_osptr.get();
    }

private:
    std::unique_ptr<std::ostream> _denovo_osptr;
};

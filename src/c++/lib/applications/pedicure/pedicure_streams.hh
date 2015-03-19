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

#include "pedicure_shared.hh"
#include "PedicureSampleSetSummary.hh"

#include "starling_common/starling_streams_base.hh"



struct pedicure_streams : public starling_streams_base
{
    typedef starling_streams_base base_t;

    pedicure_streams(
        const pedicure_options& opt,
        const pedicure_deriv_options& dopt,
        const prog_info& pinfo,
        const bam_header_t* const bam_header,
        const PedicureSampleSetSummary& ssi);

    std::ostream*
    denovo_osptr() const
    {
        return _denovo_osptr.get();
    }

    std::ostream*
    denovo_callable_osptr() const
    {
        return _denovo_callable_osptr.get();
    }

private:
    std::unique_ptr<std::ostream> _denovo_osptr;
    std::unique_ptr<std::ostream> _denovo_callable_osptr;
};

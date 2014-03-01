// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

/// \file

/// \author Chris Saunders
///

#pragma once

#include "strelka_sample_type.hh"

#include "starling_common/starling_streams_base.hh"
#include "strelka/strelka_shared.hh"



struct strelka_streams : public starling_streams_base {

    typedef starling_streams_base base_t;

    strelka_streams(const strelka_options& opt,
                    const prog_info& pinfo,
                    const bam_header_t* const bam_header);

    ~strelka_streams();

    std::ostream*
    somatic_snv_osptr() const
    {
        return _somatic_snv_osptr.get();
    }

    std::ostream*
    somatic_indel_osptr() const
    {
        return _somatic_indel_osptr.get();
    }

    std::ostream*
    somatic_callable_osptr() const
    {
        return _somatic_callable_osptr.get();
    }

private:
    std::auto_ptr<std::ostream> _somatic_snv_osptr;
    std::auto_ptr<std::ostream> _somatic_indel_osptr;
    std::auto_ptr<std::ostream> _somatic_callable_osptr;
};

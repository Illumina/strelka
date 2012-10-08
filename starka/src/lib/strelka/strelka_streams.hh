// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2012 Illumina, Inc.
//
// This software is covered by the "Illumina Genome Analyzer Software
// License Agreement" and the "Illumina Source Code License Agreement",
// and certain third party copyright/licenses, and any user of this
// source file is bound by the terms therein (see accompanying files
// Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
// Illumina_Source_Code_License_Agreement.pdf and third party
// copyright/license notices).
//
//

/// \file

/// \author Chris Saunders
///
#ifndef __STRELKA_STREAMS_HH
#define __STRELKA_STREAMS_HH

#include "strelka_sample_type.hh"

#include "starling_common/starling_streams_base.hh"
#include "strelka/strelka_shared.hh"


#define OUTPUT_VCF


struct strelka_streams : public starling_streams_base {

    typedef starling_streams_base base_t;

    strelka_streams(const strelka_options& opt,
                    const prog_info& pinfo,
                    const bam_header_t* const bam_header);

    ~strelka_streams();

    std::ostream*
    somatic_snv_osptr() const {
        return _somatic_snv_osptr.get();
    }

    std::ostream*
    somatic_indel_osptr() const {
        return _somatic_indel_osptr.get();
    }

private:
    std::auto_ptr<std::ostream> _somatic_snv_osptr;
    std::auto_ptr<std::ostream> _somatic_indel_osptr;
};


#endif

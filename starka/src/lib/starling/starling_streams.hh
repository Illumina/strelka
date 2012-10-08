// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright 2009 Illumina, Inc.
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
#ifndef __STARLING_STREAMS_HH
#define __STARLING_STREAMS_HH

#include "starling_common/starling_streams_base.hh"


struct starling_streams : public starling_streams_base {

    typedef starling_streams_base base_t;

    starling_streams(const starling_options& client_opt,
                     const prog_info& pinfo,
                     const bam_header_t* const bam_header);
};


#endif

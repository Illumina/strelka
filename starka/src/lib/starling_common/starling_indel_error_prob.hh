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
///
/// \author Chris Saunders
///


#ifndef __STARLING_INDEL_ERROR_PROB_HH
#define __STARLING_INDEL_ERROR_PROB_HH


#include "starling_common/starling_indel_report_info.hh"
#include "starling_common/starling_shared.hh"


//
// "indel_error" is the probability that the read supporting the indel case is an error
// "ref_error" is the probability that the read supporting the ref case is an error
//
void
get_indel_error_prob(const starling_options& client_opt,
                     const starling_indel_report_info& iri,
                     double& indel_error_prob,
                     double& ref_error_prob);

#endif


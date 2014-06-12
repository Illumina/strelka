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


// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/downloads/sequencing/licenses/>.
//

/// \file
///
/// \author Chris Saunders
///

#ifndef __STARLING_READ_ALIGN_SCORE_INDELS_HH
#define __STARLING_READ_ALIGN_SCORE_INDELS_HH


#include "candidate_alignment.hh"

#include "starling_common/indel_synchronizer.hh"
#include "starling_common/starling_read_segment.hh"
#include "starling_common/starling_shared.hh"

#include <set>
#include <vector>


typedef indel_buffer::iterator iiter;
typedef indel_buffer::const_iterator ciiter;

typedef std::map<indel_key,bool> indel_status_map_t;


// use the most likely alignment for each indel state for every indel
// in indel_status_map to generate data needed in indel calling:
//
void
score_indels(const starling_options& opt,
             const starling_deriv_options& dopt,
             const starling_sample_options& sample_opt,
             const read_segment& rseg,
             indel_synchronizer& isync,
             const std::set<candidate_alignment>& cal_set,
             const bool is_incomplete_search,
             const std::vector<double>& cal_set_path_lnp,
             double max_path_lnp,
             const candidate_alignment* max_cal_ptr);

#endif

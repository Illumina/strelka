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

#ifndef __ADJUST_JOINT_EPROB_HH
#define __ADJUST_JOINT_EPROB_HH

#include "blt_common/blt_shared.hh"
#include "blt_common/snp_pos_info.hh"

#include <vector>


// cache dependent probs at a fixed vexp value (cache should be
// associated with a specific blt_options object):
//
struct dependent_prob_cache {

    enum { MAX_QSCORE = 64 };

    dependent_prob_cache() : _val(MAX_QSCORE+1), _is_init(MAX_QSCORE+1,false) {}

    // client is contracted to always call this function with the same
    // vexp value used on the first call
    //
    blt_float_t
    get_dependent_val(const unsigned qscore,
                      const blt_float_t vexp);

private:
    std::vector<blt_float_t> _val;
    std::vector<bool> _is_init;
};



void
adjust_joint_eprob(const blt_options& opt,
                   dependent_prob_cache& dpc,
                   const snp_pos_info& pi,
                   std::vector<float>& dependent_eprob,
                   const bool is_depenent);

#endif

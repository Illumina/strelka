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
#ifndef __SAMPLE_RANGE_HH
#define __SAMPLE_RANGE_HH


// Interface for sample function object, sampled function should be
// able to inherit from this interface.
//
#if 0
struct sample_func_iface {
    virtual
    blt_float_t
    calc_and_store_val(const blt_float_t x) = 0;
};
#endif


// A cheap deterministic pseudo-importance sampler
//
// Object func must meet requirements described above for
// sample_func_iface
//
template <typename Func>
void
sample_uniform_range(const blt_float_t min_x,
                     const blt_float_t max_x,
                     Func& f);


#include "sample_range.hpp"

#endif

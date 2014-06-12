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

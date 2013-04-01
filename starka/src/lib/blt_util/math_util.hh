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
#ifndef __MATH_UTIL_HH
#define __MATH_UTIL_HH

#include "boost/math/special_functions/log1p.hpp"

#include <cmath>

#include <algorithm>


/// returns log(1+x), switches to special libc function when abs(x) is small
///
template <typename FloatType>
FloatType
log1p_switch(const FloatType x) {

    // better number??
    static const FloatType smallx_thresh(0.01);

    if(std::abs(x)<smallx_thresh) {
        return boost::math::log1p(x);
    } else {
        return std::log(1+x);
    }
}


/// returns equiv of log(exp(x1)+exp(x2))
///
template <typename FloatType>
FloatType
log_sum(FloatType x1, FloatType x2) {
    if(x1<x2) std::swap(x1,x2);
    return x1 + log1p_switch(std::exp(x2-x1));
}

#endif

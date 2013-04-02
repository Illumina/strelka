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

/// \author Chris Saunders
///
#ifndef __BINOMIAL_TEST_HH
#define __BINOMIAL_TEST_HH

#include <ciso646>


bool
is_reject_binomial_p_exact(const double alpha,
                           const double p,
                           const unsigned n_success,
                           const unsigned n_failure);

bool
is_reject_binomial_p_chi_sqr(const double alpha,
                             const double p,
                             const unsigned n_success,
                             const unsigned n_failure);

/// \brief test to reject a binomial with p='p' at FPR<='alpha'
///
/// chooses the testing method to strike a reasonable compromise
/// between efficiency and accuracy
///
bool
is_reject_binomial_p(const double alpha,
                     const double p,
                     const unsigned n_success,
                     const unsigned n_failure);

#endif

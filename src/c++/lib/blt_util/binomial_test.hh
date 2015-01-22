// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Starka
// Copyright (c) 2009-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

///
/// \author Chris Saunders
///

#pragma once

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

// brief test to reject a p-value from a binomial exact
// test with a target p-value.  Returns true when p-value
// is equal to or less than the threshold (i.e. null
// hypothesis is rejected) and false if the p-value is
// above threshold
// N.B. currently only evaluates a one-sided binomial
// exact test estimating the probability that a value
// less than n_success could have been generated from
// n_trials given a probability of p
bool
is_reject_binomial_pval(const double threshold,
                        const double p,
                        const unsigned n_success,
                        const unsigned n_trials);


// a one-sided binomial exact test
// Currently only estimates the probability that a value
// less than n_success could have been generated from
// n_trials given a probability of p
double
get_binomial_pval(const double p,
                  const unsigned n_success,
                  const unsigned n_trials);

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

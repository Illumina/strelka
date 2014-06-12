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

/// \file
///
/// \author Chris Saunders
///
/// routines to test independence of contingency tables
///
/// note that tables are always passed in as row-major linear arrays,
/// such that the table:
///
/// 1 2 3
/// 4 5 6
///
/// ...is represented with:
///
/// table[] =  { 1, 2, 3, 4, 5, 6 }
///

#ifndef __TABLE_TEST_HH
#define __TABLE_TEST_HH

#include <ciso646>


#ifdef HAVE_FISHER_EXACT_TEST
// optionally the workspace ws can be defined with
// get_exact_test_ws() if this is being called repeatedly
// to save on alloc/free overhead for each test
//
double
table_exact_pval(const int* table,
                 const unsigned n_row,
                 const unsigned n_col,
                 double* ws = 0);

bool
is_reject_table_exact(const double alpha,
                      const int* table,
                      const unsigned n_row,
                      const unsigned n_col,
                      double* ws = 0);

double*
get_exact_test_ws();
#endif

double
table_chi_sqr_pval(const int* table,
                   const unsigned n_row,
                   const unsigned n_col);

bool
is_reject_table_chi_sqr(const double alpha,
                        const int* table,
                        const unsigned n_row,
                        const unsigned n_col);

#if 0
/// \brief
///
/// Chooses the testing method to strike a reasonable compromise
/// between efficientcy and accuracy. Note that this test does
/// not support pre-allocing memeory for many calls to the exact
/// test, so it is not the most efficient solution for many calls!
///
bool
is_reject_table(const double alpha,
                const double* table,
                const unsigned n_row,
                const unsigned n_col);
#endif

#endif

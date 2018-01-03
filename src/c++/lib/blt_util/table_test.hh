//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2018 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

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

#pragma once


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

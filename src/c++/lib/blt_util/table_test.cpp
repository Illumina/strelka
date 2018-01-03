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

/// \file

/// \author Chris Saunders
///
#ifdef HAVE_FISHER_EXACT_TEST
extern "C" {
#include "fexact.h"
}
#endif

#include "blt_util/log.hh"
#include "blt_util/table_test.hh"

#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/chi_squared.hpp>

using boost::math::binomial;
using boost::math::cdf;
using boost::math::chi_squared;

#include <algorithm>

#ifdef HAVE_FISHER_EXACT_TEST
//#define DEBUG_EXACT
#ifdef DEBUG_EXACT
#include <iostream>
#endif

int workspace(2000000);



double*
get_exact_test_ws()
{
    double* ws = make_ws(workspace);
    if (ws==0)
    {
        log_os << "ERROR: can't allocate exact test workspace.\n";
        exit(EXIT_FAILURE);
    }
    return ws;
}



double
table_exact_pval(const int* table,
                 const unsigned n_row,
                 const unsigned n_col,
                 double* ws)
{

    static const unsigned MAX_DIM(10);

    bool is_free_ws(false);
    if (ws==0)
    {
        ws=get_exact_test_ws();
        is_free_ws=true;
    }

    assert((n_row != 0) and (n_col != 0));
    assert((n_row <= MAX_DIM) and (n_col <= MAX_DIM));

    int nr(n_row);
    int nc(n_col);

    double expect(5.);
    double percent(80.);
    double emin(1.);
    double prt(0.);
    double pval(0.);
    int mult(30);

    fexact(&nc,&nr,const_cast<int*>(table),&nc,&expect,&percent,&emin,&prt,&pval,&workspace,ws,&mult);

    if (is_free_ws) free(ws);

#ifdef DEBUG_EXACT
    std::cerr << "table:\n";
    for (unsigned i(0); i<n_row; ++i)
    {
        for (unsigned j(0); j<n_col; ++j)
        {
            if (j) std::cerr << " ";
            std::cerr << table[j+i*n_col];
        }
        std::cerr << "\n";
    }
    std::cerr "pval/alpha: " << pval << " " << alpha << "\n";
#endif

    return pval;
}



bool
is_reject_table_exact(const double alpha,
                      const int* table,
                      const unsigned n_row,
                      const unsigned n_col,
                      double* ws)
{

    return (table_exact_pval(table,n_row,n_col,ws) < alpha);
}
#endif



double
table_chi_sqr_pval(const int* table,
                   const unsigned n_row,
                   const unsigned n_col)
{

    static const unsigned MAX_DIM(10);

    assert((n_row > 1) && (n_col > 1));
    assert((n_row <= MAX_DIM) && (n_col <= MAX_DIM));

    int sum(0);
    int rsum[MAX_DIM];
    int csum[MAX_DIM];

    for (unsigned r(0); r<n_row; ++r) rsum[r] = 0;
    for (unsigned c(0); c<n_col; ++c) csum[c] = 0;

    for (unsigned r(0); r<n_row; ++r)
    {
        for (unsigned c(0); c<n_col; ++c)
        {
            const int obs(table[c+r*n_col]);
            assert(obs>=0);
            csum[c] += obs;
            rsum[r] += obs;
            sum += obs;
        }
    }

    if (sum <= 0.) return 1.;

    double xsq(0);
    for (unsigned r(0); r<n_row; ++r)
    {
        for (unsigned c(0); c<n_col; ++c)
        {
            const int obs(table[c+r*n_col]);
            const double expect((rsum[r]*csum[c])/sum);
            if (expect <= 0) continue;
            const double d(obs-expect);
            xsq += (d*d)/expect;
        }
    }

    unsigned n_df_row(n_row);
    unsigned n_df_col(n_col);
    for (unsigned r(0); r<n_row; ++r)
    {
        if (rsum[r] == 0)
        {
            n_df_row--;
            if (n_df_row<=1) return 1.;
        }
    }
    for (unsigned c(0); c<n_col; ++c)
    {
        if (csum[c] == 0)
        {
            n_df_col--;
            if (n_df_col<=1) return 1.;
        }
    }
    const unsigned df((n_df_row-1)*(n_df_col-1));
    chi_squared dist(df);
    return (1.-cdf(dist,xsq));
}



bool
is_reject_table_chi_sqr(const double alpha,
                        const int* table,
                        const unsigned n_row,
                        const unsigned n_col)
{

    return (table_chi_sqr_pval(table,n_row,n_col) < alpha);
}



#if 0
bool
is_reject_table(const double alpha,
                const int* table,
                const unsigned n_row,
                const unsigned n_col)
{

    static const unsigned exact_test_threshold(250);

    const unsigned n_trial(n_success+n_failure);

    if (n_trial > exact_test_threshold)
    {
        return is_reject_table_chi_sqr(alpha,table,n_row,n_col);
    }
    else
    {
        return is_reject_table_exact(alpha,table,n_row,n_col);
    }
}
#endif

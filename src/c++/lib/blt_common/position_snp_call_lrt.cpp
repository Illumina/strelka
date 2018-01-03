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
#include "position_loghood_minfunc.hh"
#include "blt_common/position_snp_call_lrt.hh"

#define CODEMIN_USE_BOOST
#include "minimize_1d.h"
#include "minimize_conj_direction.h"

#include "boost/math/distributions/chi_squared.hpp"

#include <cassert>

#include <iostream>



void
position_snp_call_lrt(const double alpha,
                      const snp_pos_info& pi,
                      lrt_snp_call& sc)
{
    if (pi.get_ref_base()=='N') return;

    unsigned ecount(0);
    const unsigned n_calls(pi.calls.size());
    for (unsigned i(0); i<n_calls; ++i)
    {
        assert(pi.calls[i].base_id !=BASE_ID::ANY);
        ecount++;
    }
    if (ecount==0) return;

    const double null_loghood(calc_pos_nonref_freq_loghood(pi,0.));

    // heuristic to escape early:
    static const double p_delta(0.001);
    const double delta_loghood(calc_pos_nonref_freq_loghood(pi,p_delta));
    if (null_loghood > delta_loghood) return;

    double x_nonref_freq;
    double x_loghood;

    position_nonref_freq_loghood_minfunc mf(pi);

    static const double x1(0.5);
    static const double x2(0.4);
    codemin::minimize_1d(x1,x2,mf.val(x1),mf,x_nonref_freq,x_loghood);

    x_nonref_freq = mf.arg_to_prob(x_nonref_freq);

    const double log_lrt(-2.*(x_loghood+null_loghood));

    // because null has the parameter fixed to a boundary value, the
    // asymmtotic distribution is a 50:50 mixture of csq(0) and chq(1)
    // -- the same effect as multiplying alpha of csq(1) by 2, dividing
    // the null prob by 2. (as we do below):
    boost::math::chi_squared dist(1);
    const double null_prob((1.-boost::math::cdf(dist,log_lrt))/2.);

    sc.is_snp=(null_prob<alpha);
    sc.null_loghood=null_loghood;
    sc.min_test_loghood=-x_loghood;
    sc.snp_prob=1.-null_prob;

    // if it's a snp then get additional information on non-reference
    // allele frequencies.
    //
    if (! sc.is_snp) return;

    static const double line_tol(1e-7);
    static const double start_ratio(0.05);
    static const double min_start_dist(1e-6);
    static const double end_tol(1e-7);
    static const unsigned max_iter(200);

    const unsigned ref_base_id(base_to_id(pi.get_ref_base()));

    const double ref_freq(1.-x_nonref_freq);
    const double nonref_freq((x_nonref_freq)/3.);
    for (unsigned i(0); i<N_BASE; ++i)
    {
        if (i==ref_base_id) sc.allele_freq[i] = ref_freq;
        else               sc.allele_freq[i] = nonref_freq;
    }

    static const unsigned N_BASE2(N_BASE*N_BASE);
    double conj_dir[N_BASE2];
    std::fill(conj_dir,conj_dir+N_BASE2,0.);
    for (unsigned i(0); i<N_BASE; ++i)
    {
        const double start_dist( std::max(std::fabs(sc.allele_freq[i]*start_ratio),min_start_dist) );
        conj_dir[i*(N_BASE+1)] = start_dist;
    }

    double start_tol(end_tol);
    unsigned iter;
    double x_all_loghood;
    double final_dlh;
    position_allele_distro_loghood_minfunc alm(pi);
    codemin::minimize_conj_direction(sc.allele_freq,conj_dir,alm,start_tol,end_tol,line_tol,
                                     x_all_loghood,iter,final_dlh,max_iter);
    alm.arg_to_prob(sc.allele_freq,sc.allele_freq);

    sc.min_loghood=-x_all_loghood;
}



std::ostream& operator<<(std::ostream& os,
                         const lrt_snp_call& sc)
{

    os << std::setprecision(10) << std::fixed;

    os << "P(snp): " << sc.snp_prob
       << std::setprecision(8)
       << " freq(A): " << sc.allele_freq[0] << " freq(C): "  << sc.allele_freq[1]
       << " freq(G): " << sc.allele_freq[2] << " freq(T): "  << sc.allele_freq[3];

    os.unsetf(std::ios::fixed);

    return os;
}

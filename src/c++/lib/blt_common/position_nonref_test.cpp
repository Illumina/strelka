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
#include "blt_common/position_nonref_test.hh"
#include "blt_util/math_util.hh"
#include "blt_util/prob_util.hh"

#include "position_loghood_minfunc.hh"
#include "sparse_function_1d.hh"
#include "sample_range.hh"

#define CODEMIN_USE_BOOST
#include "minimize_1d.h"
#include "minimize_conj_direction.h"

//#include "boost/math/distributions/chi_squared.hpp"

#include <cassert>

#include <iostream>



struct nonref_freq_loghood_sparse_func   //: public sample_func_iface {
{

    nonref_freq_loghood_sparse_func(const snp_pos_info& pi,
                                    sparse_function& sf)
        : _pi(pi), _sf(sf) {}

    blt_float_t
    calc_val(const blt_float_t x) const
    {
        return calc_pos_nonref_freq_loghood(_pi,x);
    }

    blt_float_t
    calc_and_store_val(const blt_float_t x)
    {
        const blt_float_t val(calc_val(x));
        _sf.insert(x,val);
        return val;
    }

private:
    const snp_pos_info& _pi;
    sparse_function& _sf;
};



void
position_nonref_test(const snp_pos_info& pi,
                     const double nonref_site_prob,
                     const double min_nonref_freq,
                     const bool is_mle_freq,
                     nonref_test_call& nrc)
{

    if (pi.get_ref_base()=='N') return;

    // add early escape test here?

    blt_float_t lhood[NRTEST::SIZE];

    lhood[NRTEST::REF] = calc_pos_nonref_freq_loghood(pi,0.);

    sparse_function sf;

#if 1
    nonref_freq_loghood_sparse_func nlf(pi,sf);
    sample_uniform_range(min_nonref_freq,1.,nlf);

#else

    // start out with a hack just to get the method rolling:
    static const blt_float_t freqs[] = { 0.01, 0.05, 0.25, 0.5, 0.75, 1. };
    static const unsigned n_freq=sizeof(freqs)/sizeof(blt_float_t);
    for (unsigned i(0); i<n_freq; ++i)
    {
        sf.insert(freqs[i],calc_pos_nonref_freq_loghood(pi,freqs[i]));
    }

#endif

    lhood[NRTEST::NONREF] = integrate_ln_sparsefunc(sf, min_nonref_freq, 1,1,1);

    //std::cerr << "WAGART: logh ref/nonef: " << lhood[0] << " " << lhood[1] << "\n";

    // TODO: ctor compute this:
    blt_float_t prior[NRTEST::SIZE];
    prior[NRTEST::REF] = log1p_switch(static_cast<double>(-nonref_site_prob));
    prior[NRTEST::NONREF] = std::log(nonref_site_prob);

    double pprob[NRTEST::SIZE];
    for (unsigned i(0); i<NRTEST::SIZE; ++i)
    {
        pprob[i] = lhood[i] + prior[i];
    }
    normalizeLogDistro(pprob, pprob + NRTEST::SIZE, nrc.max_gt);

    nrc.snp_qphred=error_prob_to_qphred(pprob[NRTEST::REF]);
    nrc.max_gt_qphred=error_prob_to_qphred(prob_comp(pprob,pprob+NRTEST::SIZE,nrc.max_gt));

    nrc.is_snp=(nrc.snp_qphred != 0);

    if (! (is_mle_freq && nrc.is_snp)) return;

#if 0
    const double null_loghood(calc_pos_nonref_freq_loghood(pi,0.));

    // heuristic to escape early:
    static const double p_delta(0.001);
    const double delta_loghood(calc_pos_nonref_freq_loghood(pi,p_delta));
    if (null_loghood > delta_loghood) return;

    double x_nonref_freq;
    double x_loghood;

    position_nonref_freq_loghood_minfunc mf(epi);

    static const double x1(0.5);
    static const double x2(0.4);
    codemin::minimize_1d(x1,x2,mf.val(x1),mf,x_nonref_freq,x_loghood);

    x_nonref_freq = mf.arg_to_prob(x_nonref_freq);

    const double log_lrt(-2.*(x_loghood+null_loghood));

    // becuase null has the parameter fixed to a boundary value, the
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
    if (not sc.is_snp) return;

    static const double line_tol(1e-7);
    static const double start_ratio(0.05);
    static const double min_start_dist(1e-6);
    static const double end_tol(1e-7);
    static const unsigned max_iter(200);

    const unsigned ref_base_id(base_to_id(pi.ref_base));

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
    position_allele_distro_loghood_minfunc alm(epi);
    codemin::minimize_conj_direction(sc.allele_freq,conj_dir,alm,start_tol,end_tol,line_tol,
                                     x_all_loghood,iter,final_dlh,max_iter);
    alm.arg_to_prob(sc.allele_freq,sc.allele_freq);

    sc.min_loghood=-x_all_loghood;
#endif
}



void
write_nonref_test(const blt_options& /*opt*/,
                  const snp_pos_info& /*pi*/,
                  const nonref_test_call& nrc,
                  std::ostream& os)
{

    os << nrc.snp_qphred
       << '\t' << NRTEST::label(static_cast<NRTEST::index_t>(nrc.max_gt))
       << '\t' << nrc.max_gt_qphred;

#if 0
    if (opt.is_print_all_poly_gt)
    {
        for (unsigned gt(0); gt<DIGT::SIZE; ++gt)
        {
#if 1
            // print GT as prob:
            os << '\t' << po.pprob[gt];
#else
            // print GT as qval:
            os << '\t' << error_prob_to_qphred(prob_comp(po.pprob,po.pprob+DIGT::SIZE,gt));
#endif
        }
    }
    const result_set& ge(dgt.genome);
    const result_set& po(dgt.poly);
#endif

#if 0
    if (nrc.is_freq)
    {
        os << std::setprecision(8) << std::fixed;
        for (unsigned i(0); i<N_BASE; ++i)
        {
            os << '\t' << nrc.allele_freq[i];
        }
        os.unsetf(std::ios::fixed);
    }
#endif
}


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
#include "blt_common/position_nonref_2allele_test.hh"
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



struct nonref_allele_freq_loghood_sparse_func   //: public sample_func_iface {
{

    nonref_allele_freq_loghood_sparse_func(const snp_pos_info& pi,
                                           const unsigned nonref_id,
                                           sparse_function& sf)
        : _pi(pi), _nonref_id(nonref_id), _sf(sf) {}

    blt_float_t
    calc_val(const blt_float_t x) const
    {
        return calc_pos_nonref_allele_freq_loghood(_pi,_nonref_id,x);
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
    const unsigned _nonref_id;
    sparse_function& _sf;
};



void
position_nonref_2allele_test(
    const snp_pos_info& pi,
    const double nonref_variant_rate,
    const double min_nonref_freq,
    const double nonref_site_error_rate,
    const double nonref_site_error_decay_freq,
    const bool /*is_always_test*/,
    nonref_test_call& nrc)
{
    static const bool is_mle_freq(false);

    if (pi.get_ref_base()=='N') return;

    // add early escape test here?

    // 1. Determine the two 'primary' alleles -- Simple test just adds
    // up qscores to determine which alleles are primary.
    //
    nrc.nonref_id=(BASE_ID::ANY);
    //unsigned nonref2_id(BASE_ID::ANY); // just ignore this value for now....
    {
        double qtot[N_BASE];
        for (unsigned i(0); i<N_BASE; ++i) qtot[i] = 0;

        const unsigned n_calls(pi.calls.size());
        for (unsigned i(0); i<n_calls; ++i)
        {
            if (pi.calls[i].base_id==BASE_ID::ANY) continue;
            qtot[pi.calls[i].base_id] += pi.calls[i].get_qscore();
        }

        // get max and max2:
        unsigned max_id=0;
        unsigned max2_id=1;
        for (unsigned b(1); b<N_BASE; ++b)
        {
            if (qtot[b] > qtot[max_id])
            {
                max2_id = max_id;
                max_id = b;
            }
            else if (qtot[b] > qtot[max2_id])
            {
                max2_id = b;
            }
        }

        const unsigned ref_id=base_to_id(pi.get_ref_base());
        if       (ref_id==max_id)
        {
            nrc.nonref_id=max2_id;

#if 0
        }
        else if (ref_id==max2_id)
        {
            nrc.nonref_id=max_id;
#endif
        }
        else
        {
            nrc.nonref_id=max_id;
            //nonref2_id=max2_id;
        }
    }

    blt_float_t lhood[NR2TEST::SIZE];

    lhood[NR2TEST::REF] = calc_pos_nonref_freq_loghood(pi,0.);

    sparse_function sf;
    nonref_allele_freq_loghood_sparse_func nlf(pi,nrc.nonref_id,sf);
    sample_uniform_range(0.,1.,nlf);
    //sample_uniform_range(min_nonref_freq,1.,nlf);

    lhood[NR2TEST::NONREF_MF] = integrate_ln_sparsefunc(sf, min_nonref_freq, 1,1,1);
    lhood[NR2TEST::NONREF_MF_NOISE] = integrate_ln_sparsefunc(sf, 0, nonref_site_error_decay_freq,1,0);

    static const blt_float_t neginf(-std::numeric_limits<blt_float_t>::infinity());
    lhood[NR2TEST::NONREF_OTHER] = neginf;

    //std::cerr << "WAGART: logh ref/nonef: " << lhood[0] << " " << lhood[1] << "\n";

    // TODO: ctor compute this:

    // this goes in here just in case someone cranks both parameters up near 1:
    //
    const double nonref_variant_rate_used = nonref_variant_rate*(1-nonref_site_error_rate);

    blt_float_t prior[NR2TEST::SIZE];
    prior[NR2TEST::REF] = log1p_switch(-(nonref_variant_rate_used+nonref_site_error_rate));
    prior[NR2TEST::NONREF_MF] = std::log(nonref_variant_rate_used/3);
    prior[NR2TEST::NONREF_MF_NOISE] = std::log(nonref_site_error_rate);
    prior[NR2TEST::NONREF_OTHER] = std::log(2*nonref_variant_rate_used/3);

    double pprob[NR2TEST::SIZE];
    for (unsigned i(0); i<NR2TEST::SIZE; ++i)
    {
        pprob[i] = lhood[i] + prior[i];
    }
    normalizeLogDistro(pprob, pprob + NR2TEST::SIZE, nrc.max_gt);

    nrc.snp_qphred=error_prob_to_qphred(pprob[NR2TEST::REF]+pprob[NR2TEST::NONREF_MF_NOISE]);
    nrc.max_gt_qphred=error_prob_to_qphred(prob_comp(pprob,pprob+NR2TEST::SIZE,nrc.max_gt));

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
write_nonref_2allele_test(const snp_pos_info& pi,
                          const nonref_test_call& nrc,
                          std::ostream& os)
{
    os << nrc.snp_qphred
       << '\t' << NR2TEST::label(static_cast<NR2TEST::index_t>(nrc.max_gt)) << "_" << id_to_base(nrc.nonref_id)
       << '\t' << nrc.max_gt_qphred;

    pi.print_known_counts(os);
    pi.print_known_qscore(os);

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

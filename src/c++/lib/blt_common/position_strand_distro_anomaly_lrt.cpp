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
#include "position_loghood_minfunc.hh"
#include "blt_common/position_strand_distro_anomaly_lrt.hh"

#include "blt_util/log.hh"
#include "blt_util/stat_util.hh"

#define CODEMIN_USE_BOOST
#include "minimize_conj_direction.h"

#include <cassert>
#include <cmath>
#include <cstdlib>

#include <iostream>


// allele map[N_BASE] is true for each base listed in allele_freq[<=N_BASE]:
//
static
void
get_max_lhood_allele_freq(const snp_pos_info& pi,
                          double* allele_freq,
                          bool* is_allele_used,
                          double& loghood)
{
    // minimization constants:
    static const double line_tol(1e-7);
    static const double start_ratio(0.05);
    static const double min_start_dist(1e-6);
    static const double end_tol(1e-7);
    static const unsigned max_iter(200);

    static const unsigned N_BASE2(N_BASE*N_BASE);
    double conj_dir[N_BASE2];

    unsigned n_allele(0);
    for (unsigned i(0); i<N_BASE; ++i)
    {
        if (is_allele_used[i]) n_allele++;
    }

    assert(n_allele);

    const unsigned n_allele2(n_allele*n_allele);

    std::fill(conj_dir,conj_dir+n_allele2,0.);
    for (unsigned i(0); i<n_allele; ++i)
    {
        const double start_dist( std::max(std::fabs(allele_freq[i]*start_ratio),min_start_dist) );
        conj_dir[i*(n_allele+1)] = start_dist;
    }

    double start_tol(end_tol);
    unsigned iter;
    double final_dlh;
    position_allele_distro_loghood_minfunc alm(pi,is_allele_used);
    codemin::minimize_conj_direction(allele_freq,conj_dir,alm,start_tol,end_tol,line_tol,
                                     loghood,iter,final_dlh,max_iter);
    alm.arg_to_prob(allele_freq,allele_freq);
    loghood=-loghood;
}



void
position_strand_distro_anomaly_lrt_expert(const snp_pos_info& pi,
                                          double& null_loghood,
                                          double& alt_loghood,
                                          unsigned& df)
{
    null_loghood=0;
    alt_loghood=0;
    df=0;

    const unsigned n_calls(pi.calls.size());
    if (n_calls == 0) return;

    bool is_allele_used[N_BASE] = {false,false,false,false};

    for (unsigned i(0); i<n_calls; ++i)
    {
        const uint8_t obs_id(pi.calls[i].base_id);
        assert(obs_id !=BASE_ID::ANY);
        is_allele_used[obs_id] = true;
    }

    unsigned n_allele(0);
    for (unsigned i(0); i<N_BASE; ++i)
    {
        if (is_allele_used[i]) n_allele++;
    }

    if (n_allele<2) return;

    snp_pos_info fstrand_pi;
    snp_pos_info rstrand_pi;

    fstrand_pi.set_ref_base(pi.get_ref_base());
    rstrand_pi.set_ref_base(pi.get_ref_base());
    for (unsigned i(0); i<n_calls; ++i)
    {
        if (pi.calls[i].is_fwd_strand)
        {
            fstrand_pi.calls.push_back(pi.calls[i]);
        }
        else
        {
            rstrand_pi.calls.push_back(pi.calls[i]);
        }
    }

    if (fstrand_pi.calls.size() == 0 || rstrand_pi.calls.size() == 0) return;

    const double allele_expect(1./static_cast<double>(n_allele));

    double joint_allele_freq[N_BASE];
    for (unsigned i(0); i<n_allele; ++i) joint_allele_freq[i] = allele_expect;

    get_max_lhood_allele_freq(pi,joint_allele_freq,is_allele_used,null_loghood);

    // refine call to see if we have any more (near) zero values after minimization:
    static const double zero_thresh(0.001);
    bool is_rerun(false);

    {
        unsigned i(0);
        while (i<n_allele)
        {
            if (joint_allele_freq[i]<zero_thresh)
            {
                is_rerun=true;
                for (unsigned j(i); (j+1)<n_allele; ++j)
                {
                    joint_allele_freq[j]=joint_allele_freq[j+1];
                }
                --n_allele;
                if (n_allele<2) return;

                unsigned allele_no(0);
                for (unsigned j(0); j<N_BASE; ++j)
                {
                    if (is_allele_used[j])
                    {
                        if (allele_no==i)
                        {
                            is_allele_used[j] = false;
                            break;
                        }
                        allele_no++;
                    }
                }
            }
            ++i;
        }
    }

    if (is_rerun)
    {
        get_max_lhood_allele_freq(pi,joint_allele_freq,is_allele_used,null_loghood);
    }

    double fstrand_loghood(0);
    double rstrand_loghood(0);
    double fstrand_allele_freq[N_BASE];
    double rstrand_allele_freq[N_BASE];
    for (unsigned i(0); i<n_allele; ++i)
    {
        fstrand_allele_freq[i] = joint_allele_freq[i];
        rstrand_allele_freq[i] = joint_allele_freq[i];
    }

    get_max_lhood_allele_freq(fstrand_pi,fstrand_allele_freq,is_allele_used,fstrand_loghood);
    get_max_lhood_allele_freq(rstrand_pi,rstrand_allele_freq,is_allele_used,rstrand_loghood);

    alt_loghood=fstrand_loghood+rstrand_loghood;
    df=n_allele-1;
}



bool
position_strand_distro_anomaly_lrt(const double alpha,
                                   const snp_pos_info& pi)
{
    double null_loghood(0);
    double alt_loghood(0);
    unsigned df(0);
    position_strand_distro_anomaly_lrt_expert(pi,null_loghood,alt_loghood,df);
    if (df == 0) return false;

    return is_lrt_reject_null(null_loghood,alt_loghood,df,alpha);
}

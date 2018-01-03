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

#include "blt_common/adjust_joint_eprob.hh"

#include <cmath>

#include <algorithm>
#include <iostream>
#include <vector>


typedef unsigned icall_t;
typedef std::vector<icall_t> icalls_t;


namespace
{

struct sort_icall_by_eprob
{
    explicit
    sort_icall_by_eprob(const snp_pos_info& pi) : _pi(pi) {}

    bool
    operator()(const icall_t& a, const icall_t& b) const
    {
        return (_pi.calls[a].get_qscore() > _pi.calls[b].get_qscore());
    }

    const snp_pos_info& _pi;
};

}


static
blt_float_t
get_dependent_eprob(const unsigned qscore,
                    const blt_float_t vexp)
{
    static const blt_float_t dep_converge_prob(0.75);

    const blt_float_t eprob(qphred_to_error_prob(static_cast<int>(qscore)));
    const blt_float_t val(std::pow(eprob,vexp));
    const blt_float_t frac((1-val)/(1-eprob));
    return std::max(eprob,frac*val+(1-frac)*dep_converge_prob);
}



blt_float_t
dependent_prob_cache::
get_dependent_val(const unsigned qscore,
                  const blt_float_t vexp)
{
    assert(qscore <= MAX_QSCORE);
    if (! _is_init[qscore])
    {
        _val[qscore] = get_dependent_eprob(qscore,vexp);
        _is_init[qscore] = true;
    }

    return _val[qscore];
}



//#define DEBUG_ADJUST

static
void
adjust_icalls_eprob(const blt_options& opt,
                    dependent_prob_cache& dpc,
                    icalls_t& ic,
                    const snp_pos_info& pi,
                    std::vector<float>& dependent_eprob)
{
    const unsigned ic_size(ic.size());

#ifdef DEBUG_ADJUST
    for (unsigned i(0); i<ic_size; ++i)
    {
        base_call& bi(pi.calls[ic[i]]);
        std::cerr << "BEFORE: " << i << " " << bi.is_neighbor_mismatch << " " << dependent_eprob[ic[i]] << "\n";
    }
#endif

    static const bool is_use_vexp_frac(true);

    // produce weighted fraction of reads with a neighboring mismatch:
    blt_float_t vexp_frac;
    if (is_use_vexp_frac)
    {
        static const blt_float_t lnran(std::log(0.75));
        blt_float_t num(0);
        blt_float_t den(0);
        for (unsigned i(0); i<ic_size; ++i)
        {
            const base_call& bi(pi.calls[ic[i]]);
//            const blt_float_t eprob(dependent_eprob[ic[i]]);
            const blt_float_t weight(lnran-bi.ln_error_prob());
            den += weight;
            if (bi.is_neighbor_mismatch)
            {
                num += weight;
            }
        }
        blt_float_t mismatch_frac(0);
        if (ic_size && (den>0.)) mismatch_frac=(num/den);
        vexp_frac=(1-mismatch_frac)*opt.bsnp_ssd_no_mismatch+mismatch_frac*opt.bsnp_ssd_one_mismatch;
    }
    else
    {
        vexp_frac=opt.bsnp_ssd_no_mismatch;
    }

    const bool is_limit_vexp(opt.is_min_vexp);
    const blt_float_t min_vexp(opt.min_vexp);

    // used cached dependent probs once we reach the min_vexp level:
    bool is_min_vexp(false);

    std::sort(ic.begin(),ic.end(),sort_icall_by_eprob(pi));
    blt_float_t vexp(1.);
    for (unsigned i(0); i<ic_size; ++i)
    {
        const base_call& bi(pi.calls[ic[i]]);
        if (! is_min_vexp)
        {
            dependent_eprob[ic[i]] = static_cast<float>(get_dependent_eprob(bi.get_qscore(),vexp));

            blt_float_t next_vexp(vexp);
            if (is_use_vexp_frac)
            {
                next_vexp *= (1-vexp_frac);
            }
            else
            {
                if (bi.is_neighbor_mismatch)
                {
                    next_vexp *= (1-opt.bsnp_ssd_one_mismatch);
                }
                else
                {
                    next_vexp *= (1-opt.bsnp_ssd_no_mismatch);
                }
            }
            if (is_limit_vexp)
            {
                is_min_vexp=(next_vexp<=min_vexp);
                vexp = std::max(min_vexp,next_vexp);
            }
            else
            {
                vexp = next_vexp;
            }
        }
        else
        {
            // cached version:
            dependent_eprob[ic[i]] = static_cast<float>(dpc.get_dependent_val(bi.get_qscore(),vexp));
        }
    }

#ifdef DEBUG_ADJUST
    for (unsigned i(0); i<ic_size; ++i)
    {
        base_call& bi(pi.calls[ic[i]]);
        std::cerr << "AFTER: " << i << " " << bi.is_neighbor_mismatch << " " << dependent_eprob[ic[i]] << "\n";
    }
#endif
}



// downweight qvals at a site based on dependency effects:
//
void
adjust_joint_eprob(const blt_options& opt,
                   dependent_prob_cache& dpc,
                   const snp_pos_info& pi,
                   std::vector<float>& dependent_eprob)
{
    const unsigned n_calls(pi.calls.size());
    dependent_eprob.clear();
    for (unsigned i(0); i<n_calls; ++i)
    {
        dependent_eprob.push_back(static_cast<float>(pi.calls[i].error_prob()));
    }

    if (! opt.is_dependent_eprob()) return;

    // split icalls into fwd and reverse strand and allele types:
    //
    static const unsigned group_size(8); // (is_fwd*base_id)
    icalls_t icalls[group_size];
    for (unsigned i(0); i<n_calls; ++i)
    {
        // exclude q2's and filtered bases:
        const base_call& b(pi.calls[i]);
#ifdef NOREFFILTER
        if (b.is_call_filter && (b.base_id != pi.ref_base_id)) continue;
#else
        if (b.is_call_filter) continue;
#endif

        //if(b.is_call_filter or (b.error_prob>=0.5)) continue;
        if (b.get_qscore()<3) continue;

        const unsigned group_index((b.is_fwd_strand)+(2*b.base_id));
        icalls[group_index].push_back(i);
    }

    // process each isize array:
    for (unsigned i(0); i<group_size; ++i)
    {
        adjust_icalls_eprob(opt,dpc,icalls[i],pi,dependent_eprob);
    }
}



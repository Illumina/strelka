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

#include "blt_util/log.hh"

#include <cassert>
#include <cmath>
#include <cstdlib>

#include <iostream>


const double one_third(1./3.);



/////////////////////////////////////////////////////////////////////////////////////////////
//
// simple non-ref vs ref frequency case:
//


///
/// likelihood function is P(obs|model)
/// where the model parameter is the reference allele frequency
///
/// this becomes:
/// sum{ s in seqs } { P(obs|s) * P(s|model) }
///
/// a direct expansion would look like:
/// L = P(obs=TC|true=AA) * P(true=AA|freq=.9)+
///     P(obs=TC|true=AC) * P(true=AC|freq=.9)+...
///
/// L = P(e1)/3*P(e2)/3*.1/3*.1/3 +
///     P(e1)/3*(1-P(e2))*.1/3*.9 +...
///
/// ...where P(e1) is the error prob for basecall 1 and P(e2) is the error
/// prob for basecall 2.
///
/// The above sequence can be factored into:
/// prod { pos in column_size } { P(A|pos,model)+P(C|pos,model)+... }
///
///
///
/// The terms used below for individual basecalls are found as follows.
/// f is nonref frequency, e is basecall error prob:
///
/// Ref case:
/// Ref base prob is: P(ref) = (1-e)(1-f)
/// Non-ref base prob is: P(nonref) = (e/3)(f/3)
///
/// P(ref) + 3*P(nonref) = (1-e)(1-f)+(e)(f)/3
///
///
/// NonRef case:
/// P(ref) = (e/3)(1-f)
/// P(nonref_match) = (1-e)(f/3)
/// P(nonref_mismatch) = (e/3)(f/3)
///
/// P(ref) + P(nonref_match) + 2*P(nonref_mismatch) = (e/3)(1-f) + ((1-e)+2(e/3))(f/3)
///
blt_float_t
calc_pos_nonref_freq_loghood(const snp_pos_info& pi,
                             const blt_float_t nonref_freq)
{

    if ((nonref_freq < 0.) || (nonref_freq > 1.))
    {
        log_os << "ERROR: Invalid probability value: " << nonref_freq << "\n";
        exit(EXIT_FAILURE);
    }

    // use a simple function-local cache:
    static const unsigned max_qscore(qphred_cache::MAX_QSCORE);
    std::vector<bool> is_cache_ref(max_qscore+1,false);
    std::vector<bool> is_cache_nonref(max_qscore+1,false);

    std::vector<bool>* is_cache[2] = { &is_cache_ref, &is_cache_nonref };
    blt_float_t cache[2][max_qscore+1];


    const unsigned n_calls(pi.calls.size());

    blt_float_t loghood(0);

    const uint8_t ref_id(base_to_id(pi.get_ref_base()));
    const blt_float_t ref_freq(1-nonref_freq);

    for (unsigned i(0); i<n_calls; ++i)
    {
        const unsigned qscore(pi.calls[i].get_qscore());
        const bool is_ref_obs(ref_id == pi.calls[i].base_id);

        assert(qscore<=max_qscore);

        if (! is_cache[is_ref_obs]->operator[](qscore) )
        {
            const double eprob(qphred_to_error_prob(qscore));
            const double ceprob(1-eprob);

            // each term is a sum of the probability of each of the four bases
            if (is_ref_obs)
            {
                cache[is_ref_obs][qscore] = std::log((ceprob)*(ref_freq)+eprob*(nonref_freq*one_third));
            }
            else
            {
                cache[is_ref_obs][qscore] = std::log((eprob*one_third)*(ref_freq)+(ceprob+2*eprob*one_third)*(nonref_freq*one_third));
            }

            is_cache[is_ref_obs]->operator[](qscore) = true;
        }

        loghood += cache[is_ref_obs][qscore];
    }

    return loghood;
}



double
position_nonref_freq_loghood_minfunc::
val(const double nonref_freq) const
{
    return -calc_pos_nonref_freq_loghood(_pi,arg_to_prob(nonref_freq));
}



// a sawtooth function to keep the argument in range, yet continuous:
double
position_nonref_freq_loghood_minfunc::
arg_to_prob(const double arg)
{
    double prob(std::fabs(arg));
    const double pf(std::floor(prob));
    prob-=pf;
    if (1==static_cast<int>(pf)%2)
    {
        prob=1.-prob;
    }
    return prob;
}








/////////////////////////////////////////////////////////////////////////////////////////////
//
// simple non-ref vs specific ref frequency case:
//


///
/// The terms used below for individual basecalls are found as follows.
/// f is nonref frequency, e is basecall error prob:
///
/// Ref case:
/// Ref base prob is: P(ref) = (1-e)(1-f)
/// Id'd Non-ref base prob is: P(nonref_id) = (e/3)(f)
/// Other Non-ref base prob is: P(nonref_other) = (e/3)(0)
///
/// P(ref) + P(id) + 2*P(other) = (1-e)(1-f)+(e)(f)/3
///
///
/// NonRef_id case:
/// P(ref) = (e/3)(1-f)
/// P(nonref_id) = (1-e)(f)
/// P(nonref_mismatch) = (e/3)(0)
///
/// P(ref) + P(id) + 2*P(other) = (e/3)(1-f) + (1-e)(f)
///
/// NonRef_other case:
/// P(ref) = (e/3)(1-f)
/// P(nonref_id) = (e/3)(f)
/// P(nonref_other_match) = (1-e)(0)
/// P(nonref_other_mismatch) = (e/3)(0)
///
/// P(ref) + P(id) + P(other_match) + P(other_mismatch) = (e/3)(1-2f/3)
///
blt_float_t
calc_pos_nonref_allele_freq_loghood(const snp_pos_info& pi,
                                    const unsigned nonref_id,
                                    const blt_float_t nonref_freq)
{

    if ((nonref_freq < 0.) || (nonref_freq > 1.))
    {
        log_os << "ERROR:: invalid probability value: " << nonref_freq << "\n";
        exit(EXIT_FAILURE);
    }

    // use a simple function-local cache:
    static const unsigned max_qscore(qphred_cache::MAX_QSCORE);
    std::vector<bool> is_cache_ref(max_qscore+1,false);
    std::vector<bool> is_cache_nonref(max_qscore+1,false);
    std::vector<bool> is_cache_other(max_qscore+1,false);

    std::vector<bool>* is_cache[3] = { &is_cache_ref, &is_cache_nonref, &is_cache_other };
    blt_float_t cache[3][max_qscore+1];


    const unsigned n_calls(pi.calls.size());

    blt_float_t loghood(0.);

    const uint8_t ref_id(base_to_id(pi.get_ref_base()));
    const blt_float_t ref_freq(1.-nonref_freq);

    for (unsigned i(0); i<n_calls; ++i)
    {
        const base_call& bc(pi.calls[i]);
        const unsigned qscore(bc.get_qscore());

        const uint8_t base_id_state((ref_id == bc.base_id) ? 0 : ( (nonref_id == bc.base_id) ? 1 : 2));

        assert(qscore<=max_qscore);

        if (! is_cache[base_id_state]->operator[](qscore) )
        {
            const double eprob(qphred_to_error_prob(qscore));
            const double ceprob(1-eprob);

            if (0==base_id_state)
            {
                cache[base_id_state][qscore] = std::log((ceprob*ref_freq) + (eprob*nonref_freq*one_third));
            }
            else if (1==base_id_state)
            {
                cache[base_id_state][qscore] = std::log((eprob*one_third*ref_freq) + (ceprob*nonref_freq));
            }
            else if (2==base_id_state)
            {
                cache[base_id_state][qscore] = std::log((eprob*one_third)*(1. - 2*nonref_freq*one_third));
            }
            else
            {
                assert(0);
            }

            is_cache[base_id_state]->operator[](qscore) = true;
        }

        loghood += cache[base_id_state][qscore];
    }

    return loghood;
}








/// same as above but allow separate frequencies for all nucleotides:
///
static
double
calc_pos_allele_distro_loghood(const snp_pos_info& pi,
                               const double* allele_distro,
                               const unsigned n_allele,
                               const unsigned* allele_map)
{

    // minimization parameters should already be normalized:
    for (unsigned i(0); i<n_allele; ++i)
    {
        if ((allele_distro[i] < 0.) || (allele_distro[i] > 1.))
        {
            log_os << "ERROR: Invalid probability value: " << allele_distro[i] << "\n";
            exit(EXIT_FAILURE);
        }
    }

    const unsigned n_calls(pi.calls.size());

    double loghood(0.);

    for (unsigned i(0); i<n_calls; ++i)
    {
        const double eprob(pi.calls[i].error_prob());
        const double ceprob(1-eprob);

        const unsigned idx(pi.calls[i].base_id);
        double site_prob(0.);
        for (unsigned n(0); n<n_allele; ++n)
        {
            const double nprob( (allele_map[n] == idx) ? (ceprob) : (eprob*one_third) );
            site_prob += nprob*allele_distro[n];
        }
        loghood += std::log(site_prob);
    }

    return loghood;
}


double
position_allele_distro_loghood_minfunc::
val(const double* allele_distro_in)
{
    double allele_distro[N_BASE];
    arg_to_prob(allele_distro_in,allele_distro);
    return -calc_pos_allele_distro_loghood(_pi,allele_distro,_n_allele,_allele_map);
}



void
position_allele_distro_loghood_minfunc::
arg_to_prob(const double* allele_distro_in,
            double* allele_distro)
{

    double sum(0.);
    for (unsigned i(0); i<_n_allele; ++i)
    {
        allele_distro[i] = std::fabs(allele_distro_in[i]);
        sum += allele_distro[i];
    }

    assert(sum>=0.);
    if (sum>0.)
    {
        const double scale(1./sum);
        for (unsigned i(0); i<_n_allele; ++i) allele_distro[i] *= scale;
    }
    else
    {
        static const double allele_expect(1./static_cast<double>(_n_allele));
        for (unsigned i(0); i<_n_allele; ++i) allele_distro[i] = allele_expect;
    }
}


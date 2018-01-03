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

#pragma once

#include "blt_common/snp_pos_info.hh"
#include "blt_util/seq_util.hh"

#include <minfunc_interface.h>




/// assumes pi has already had bad cases filtered out:
///
blt_float_t
calc_pos_nonref_freq_loghood(const snp_pos_info& pi,
                             const blt_float_t nonref_freq);


// provides lhood of the reference vs. non-reference frequency:
//
struct position_nonref_freq_loghood_minfunc : public codemin::minfunc_1d_interface<double>
{

    /// assumes pi has already had bad cases filtered out:
    explicit
    position_nonref_freq_loghood_minfunc(const snp_pos_info& pi) : _pi(pi) {}

    virtual double val(const double nonref_freq) const;

    // a sawtooth function to keep the argument in range, yet continous:
    static
    double arg_to_prob(const double arg);

private:
    const snp_pos_info& _pi;
};






/// assumes pi has already had bad cases filtered out:
///
blt_float_t
calc_pos_nonref_allele_freq_loghood(const snp_pos_info& pi,
                                    const unsigned nonref_id,
                                    const blt_float_t nonref_freq);



#if 0
// provides lhood of the reference frequency vs. another allele
// frequency (or frequency of two nonref alleles):
//
struct position_nonref_allele_freq_loghood_minfunc : public codemin::minfunc_1d_interface<double>
{

    /// assumes pi has already had bad cases filtered out:
    explicit
    position_nonref_allele_freq_loghood_minfunc(const snp_pos_info& pi,
                                                const unsigned nonref_id)
        : _pi(pi)
        , _nonref_id(nonref_id)
        , _nonref2_id(nonref2_id)
    {}

    virtual double val(const double nonref_freq) const;

    // a sawtooth function to keep the argument in range, yet continous:
    static
    double arg_to_prob(const double arg);

private:
    const snp_pos_info& _pi;
    const unsigned _nonref_id;
    const unsigned _nonref2_id;
};
#endif





// provides lhood given all four allele frequencies:
//
struct position_allele_distro_loghood_minfunc : public codemin::minfunc_interface<double>
{

    // assumes pi has had bad cases filtered out already:
    // is_allele_used is true for each base submitted to val in the
    // frequency array, all other frequencies are locked at zero
    explicit
    position_allele_distro_loghood_minfunc(const snp_pos_info& pi,
                                           const bool* is_allele_used = 0) : _pi(pi), _n_allele(0)
    {
        for (unsigned i(0); i<N_BASE; ++i)
        {
            if ((! is_allele_used) || is_allele_used[i])
            {
                _allele_map[_n_allele] = i;
                _n_allele++;
            }
        }
    }

    virtual unsigned dim() const
    {
        return _n_allele;
    }

    virtual double val(const double* allele_distro_in);

    /// safe for (allele_distro_in == allele_distro) case
    void arg_to_prob(const double* allele_distro_in,
                     double* allele_distro);

private:
    const snp_pos_info& _pi;
    unsigned _n_allele;
    unsigned _allele_map[N_BASE];
};

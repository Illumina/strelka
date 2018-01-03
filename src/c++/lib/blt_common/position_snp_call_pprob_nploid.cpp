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

#include "blt_common/position_snp_call_pprob_nploid.hh"

#include "blt_util/log.hh"

#include <cassert>
#include <cmath>
#include <cstdlib>

#include <iostream>
#include <iomanip>



void
nploid_write(const nploid_info& ninfo,
             const nploid_genotype& ngt,
             std::ostream& os)
{
    os << std::setprecision(10) << std::fixed;

    os << "P(snp): " << (1.-ngt.pprob[ngt.ref_gt])
       << " max_gtype: " << ninfo.label(ngt.max_gt) << " P(max_gtype): "  << ngt.pprob[ngt.max_gt]
       << " max2_gtype: " << ninfo.label(ngt.max2_gt) << " P(max2_gtype): "  << ngt.pprob[ngt.max2_gt];

    os.unsetf(std::ios::fixed);
}



///
/// The likelihood of an observed 'column' of observations for any
/// genotype: P(obs=ACT|genotype=AT) is a sum over all possible
/// columns:
///
/// sum { c in {A,C,G,T}**3 } { P(obs=ACT|col=c) * P(col=c|genotype=AT) }
///
/// Enumerating all {A,C,G,T}**col_size possible column values is
/// inefficent, so we only consider those columns for which
/// P(col|genotype) is nonzero.  For homozygous genotypes this always
/// leaves a single column; for heterozygous genotypes this is a sum
/// over all valid columns:
///
/// P(obs=ACT|genotype=AT)=
/// sum { c in {A,T}**3 } { P(obs=ACT|col=c) * P(col=c|genotype=AT }
///
/// that is...
///
/// P(obs|AAA)*P(AAA|genotype)+
/// P(obs|AAT)*P(AAT|genotype)+...etc
///
/// the first term is a product based on error probabilities:
/// P(obs=ACT|AAA) = (1-P(e1))*(P(e2)/3)*(P(e3)/3)
///
/// the second term is a product of genotype base frequencies:
/// P(AAA|genotype=AT) = .5**3
///
/// The likelihood calculated in the function below is an equivelent
/// factorization of the calculation described above.
///
void
position_snp_call_pprob_nploid(const double snp_prob,
                               const snp_pos_info& pi,
                               const nploid_info& ninfo,
                               nploid_genotype& ngt)
{
    if (pi.get_ref_base()=='N') return;

    const unsigned n_calls(pi.calls.size());
    const unsigned ref_id(base_to_id(pi.get_ref_base()));

    // check that a non-reference call meeting quality criteria even exists:
    bool is_test(false);
    for (unsigned i(0); i<n_calls; ++i)
    {
        const uint8_t obs_id(pi.calls[i].base_id);
        assert(obs_id!=BASE_ID::ANY);
        if (ref_id!=obs_id)
        {
            is_test=true;
            break;
        }
    }

    if (! is_test) return;

    ngt.ref_gt=ninfo.get_ref_gtype(pi.get_ref_base());

    const unsigned n_gt(ninfo.gtype_size());

    // get likelihood of each genotype
    std::vector<double> lhood(n_gt,0.);

    static const double one_third(1./3.);

    const double freq_chunk(ninfo.expect_freq_chunk());
    const unsigned n_freq(ninfo.expect_freq_level_size());
    std::vector<double> ln_obs_prob_cache(n_freq);

    for (unsigned i(0); i<n_calls; ++i)
    {
        const double eprob(pi.calls[i].error_prob());

        for (unsigned j(0); j<n_freq; ++j)
        {
            const double obs_expect(j*freq_chunk);
            const double obs_prob((obs_expect)*(1.-eprob)+(1.-obs_expect)*(eprob*one_third));
            ln_obs_prob_cache[j] = std::log(obs_prob);
        }

        const uint8_t obs_id(pi.calls[i].base_id);
        for (unsigned gt(0); gt<n_gt; ++gt)
        {
            lhood[gt] += ln_obs_prob_cache[ninfo.expect_freq_level(gt,obs_id)];
        }
    }


    std::vector<double> prior(n_gt,0.);

    const double nonref_prob(snp_prob/static_cast<double>(n_gt-1));
    for (unsigned gt(0); gt<n_gt; ++gt)
    {
        if (gt==ngt.ref_gt)
        {
            prior[gt] = 1.-snp_prob;
        }
        else
        {
            prior[gt] = nonref_prob;
        }
    }

    // mult by prior distro to get unnormalized pprob:
    for (unsigned gt(0); gt<n_gt; ++gt)
    {
        ngt.pprob[gt] = lhood[gt] + std::log(prior[gt]);
    }

    // scale and exp pprob values:
    ngt.max_gt=0;
    ngt.max2_gt=1;
    double max(ngt.pprob[ngt.max_gt]);
    double max2(ngt.pprob[ngt.max2_gt]);
    for (unsigned gt(1); gt<n_gt; ++gt)
    {
        if (ngt.pprob[gt] > max)
        {
            max2 = max;
            max = ngt.pprob[gt];
            ngt.max2_gt = ngt.max_gt;
            ngt.max_gt = gt;
        }
        else if (ngt.pprob[gt] > max2)
        {
            max2 = ngt.pprob[gt];
            ngt.max2_gt = gt;
        }
    }

    // \todo debatable call criteria:
    ngt.is_snp=(ngt.max_gt != ngt.ref_gt);

    double sum(0.);
    for (unsigned gt(0); gt<n_gt; ++gt)
    {
        ngt.pprob[gt] = std::exp(ngt.pprob[gt]-max);
        sum += ngt.pprob[gt];
    }

    // normalize:
    sum = 1./sum;
    for (unsigned gt(0); gt<n_gt; ++gt)
    {
        ngt.pprob[gt] *= sum;
    }
}

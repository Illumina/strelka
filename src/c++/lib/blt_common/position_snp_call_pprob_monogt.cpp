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
#include "blt_common/position_snp_call_pprob_monogt.hh"

#include "blt_util/log.hh"
#include "blt_util/seq_util.hh"

#include <cassert>
#include <cmath>
#include <cstdlib>

#include <iostream>
#include <iomanip>



void
position_snp_call_pprob_monogt(const double theta,
                               const snp_pos_info& pi,
                               monoploid_genotype& mgt)
{
    if (pi.get_ref_base()=='N') return;

    const unsigned n_calls(pi.calls.size());
    mgt.ref_gt=base_to_id(pi.get_ref_base());

    // check that a non-reference call meeting quality criteria even exists:
    bool is_test(false);
    for (unsigned i(0); i<n_calls; ++i)
    {
        const uint8_t obs_id(pi.calls[i].base_id);
        assert(obs_id!=BASE_ID::ANY);
        if (mgt.ref_gt!=obs_id)
        {
            is_test=true;
            break;
        }
    }

    if (! is_test) return;

    // get likelihood of each genotype
    double lhood[MONOGT::SIZE];
    for (unsigned gt(0); gt<MONOGT::SIZE; ++gt) lhood[gt] = 0.;

    static const double one_third(1./3.);
    static const double log_one_third(std::log(one_third));

    for (unsigned i(0); i<n_calls; ++i)
    {
        const uint8_t obs_id(pi.calls[i].base_id);
        const double eprob(pi.calls[i].error_prob());

        const double val0(std::log(eprob)+log_one_third);
        const double val1(std::log(1.-eprob));

        for (unsigned gt(0); gt<MONOGT::SIZE; ++gt)
        {
            if (obs_id != gt)
            {
                lhood[gt] += val0;
            }
            else
            {
                lhood[gt] += val1;
            }
        }
    }

    //
    double prior[MONOGT::SIZE];
    for (unsigned gt(0); gt<MONOGT::SIZE; ++gt) prior[gt] = 0.;

    for (unsigned gt(0); gt<MONOGT::SIZE; ++gt)
    {
        if (gt==mgt.ref_gt)
        {
            prior[gt] = 1.-theta;
            assert(prior[gt]>0.);
        }
        else
        {
            prior[gt]=theta*one_third;
        }
    }

    // mult by prior distro to get unnormalized pprob:
    for (unsigned gt(0); gt<MONOGT::SIZE; ++gt)
    {
        mgt.pprob[gt] = lhood[gt] + std::log(prior[gt]);
    }

    // scale and exp pprob values:
    mgt.max_gt=0;
    mgt.max2_gt=1;
    double max(mgt.pprob[mgt.max_gt]);
    double max2(mgt.pprob[mgt.max2_gt]);
    for (unsigned gt(1); gt<MONOGT::SIZE; ++gt)
    {
        if (mgt.pprob[gt] > max)
        {
            max2 = max;
            max = mgt.pprob[gt];
            mgt.max2_gt = mgt.max_gt;
            mgt.max_gt = gt;
        }
        else if (mgt.pprob[gt] > max2)
        {
            max2 = mgt.pprob[gt];
            mgt.max2_gt = gt;
        }
    }

    // \todo debatable call criteria:
    mgt.is_snp=(mgt.max_gt != mgt.ref_gt);

    double sum(0.);
    for (unsigned gt(0); gt<MONOGT::SIZE; ++gt)
    {
        mgt.pprob[gt] = std::exp(mgt.pprob[gt]-max);
        sum += mgt.pprob[gt];
    }

    // normalize:
    sum = 1./sum;
    for (unsigned gt(0); gt<MONOGT::SIZE; ++gt)
    {
        mgt.pprob[gt] *= sum;
    }
}



std::ostream& operator<<(std::ostream& os,
                         monoploid_genotype& mgt)
{

    os << std::setprecision(10) << std::fixed;

    os << "P(snp): " << (1.-mgt.pprob[mgt.ref_gt])
       << " max_gtype: " << MONOGT::label(mgt.max_gt) << " P(max_gtype): " << mgt.pprob[mgt.max_gt]
       << " max2_gtype: " << MONOGT::label(mgt.max2_gt) << " P(max2_gtype): " << mgt.pprob[mgt.max2_gt];

    os.unsetf(std::ios::fixed);

    return os;
}

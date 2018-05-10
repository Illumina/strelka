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

#include "blt_common/snp_pos_info.hh"
#include "blt_util/io_util.hh"

#include <iomanip>
#include <iostream>
#include <type_traits>

#ifndef BC_DEBUG
static_assert(sizeof(base_call)==2, "Unexpected base_call object size");
#endif


const unsigned base_call::qscore_max = ((1 << base_call::qscore_bits) - 1);


std::ostream&
operator<<(std::ostream& os,
           const base_call& bc)
{
    const char strand(bc.is_fwd_strand ? 'F' : 'R');
    {
        const StreamScoper ss(os);
        os << "base: " << id_to_base(bc.base_id)
           << " P(error) " << std::setw(14) << std::setprecision(8) << std::left << bc.error_prob()
           << " strand: " << strand
#ifdef BC_DEBUG
           << " read_pos: " << (bc.read_pos+1)
           << " read_size: " << bc.read_size
#endif
           << " is_call_filter: " << bc.is_call_filter
           << " is_tscf?: " << bc.is_tier_specific_call_filter;
    }
    return os;
}



std::ostream&
operator<<(std::ostream& os,
           const snp_pos_info& pci)
{
    os << "ref: " << pci.get_ref_base();
    const unsigned bs(pci.calls.size());
    for (unsigned i(0); i<bs; ++i) os << "\nt1: " << pci.calls[i];
    const unsigned bs2(pci.tier2_calls.size());
    for (unsigned i(0); i<bs2; ++i) os << "\nt2: " << pci.tier2_calls[i];
    return os;
}



double
snp_pos_info::
get_read_pos_ranksum() const
{
    //	cout << read_pos_ranksum << endl;
    return readPositionRankSum.get_z_stat();

}

double
snp_pos_info::
get_mq_ranksum() const
{
    return mq_ranksum.get_z_stat();
}

double
snp_pos_info::
get_baseq_ranksum() const
{
//	cout << baseq_ranksum << endl;
    return baseq_ranksum.get_z_stat();
}

double
snp_pos_info::
get_raw_pos() const
{
//  cout << baseq_ranksum << endl;
    return readPositionRankSum.getExpectedCategory2Value();
}


double
snp_pos_info::
get_raw_baseQ() const
{
//  cout << baseq_ranksum << endl;
    return baseq_ranksum.getExpectedCategory2Value();
}


void
snp_pos_info::
print_known_counts(std::ostream& os,
                   const int min_qscore) const
{
    std::array<unsigned,N_BASE> base_count;
    getBasecallCounts(base_count, min_qscore);

    for (const unsigned bc : base_count)
    {
        os << '\t' << bc;
    }
}


void
snp_pos_info::
print_known_qscore(std::ostream& os,
                   const int min_qscore) const
{
    double qscore_tot[N_BASE];
    for (unsigned i(0); i<N_BASE; ++i) qscore_tot[i] = 0;
    unsigned qscore_count[N_BASE];
    for (unsigned i(0); i<N_BASE; ++i) qscore_count[i] = 0;

    for (const auto& call : calls)
    {
        if (call.base_id==BASE_ID::ANY) continue;
        if (call.get_qscore()<min_qscore) continue;
        qscore_tot[call.base_id] += call.get_qscore();
        os << call << '\n';
        qscore_count[call.base_id]++;
    }

    const StreamScoper ss(os);
    os << std::setprecision(2) << std::fixed;
    for (unsigned b(0); b<N_BASE; ++b)
    {
        if (qscore_count[b])
        {
            os << '\t' << qscore_tot[b]/qscore_count[b];
        }
        else
        {
            os << "\tNA";
        }
    }
}

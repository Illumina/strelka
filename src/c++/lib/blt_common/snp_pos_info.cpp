// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Starka
// Copyright (c) 2009-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

/// \file

/// \author Chris Saunders
///

#include "blt_common/snp_pos_info.hh"
#include <iomanip>
#include <iostream>



std::ostream&
operator<<(std::ostream& os,
           const base_call& bc)
{
    const char strand(bc.is_fwd_strand ? 'F' : 'R');
    os << "base: " << id_to_base(bc.base_id)
       << " P(error) " << std::setw(14) << std::setprecision(8) << std::left << bc.error_prob()
       << " strand: " << strand
#ifdef BC_DEBUG
       << " read_pos: " << (bc.read_pos+1)
       << " read_size: " << bc.read_size
#endif
       << " is_call_filter: " << bc.is_call_filter
       << " is_tscf?: " << bc.is_tier_specific_call_filter;
    return os;
}



std::ostream&
operator<<(std::ostream& os,
           const snp_pos_info& pci)
{
    os << "ref: " << pci.ref_base;
    const unsigned bs(pci.calls.size());
    for (unsigned i(0); i<bs; ++i) os << "\nt1: " << pci.calls[i];
    const unsigned bs2(pci.tier2_calls.size());
    for (unsigned i(0); i<bs2; ++i) os << "\nt2: " << pci.tier2_calls[i];
    return os;
}


double
snp_pos_info::
get_rms_mq()
{
    if (n_mapq==0)
        return 0.0;
    return std::sqrt(cumm_mapq/n_mapq);
}


double
snp_pos_info::
get_read_pos_ranksum()
{
    //	cout << read_pos_ranksum << endl;
    return read_pos_ranksum.get_u_stat();

}

double
snp_pos_info::
get_mq_ranksum()
{
//	cout << mq_ranksum << endl;
    return mq_ranksum.get_u_stat();
}

double
snp_pos_info::
get_baseq_ranksum()
{
//	cout << baseq_ranksum << endl;
    return baseq_ranksum.get_u_stat();
}

double
snp_pos_info::
get_raw_pos()
{
//  cout << baseq_ranksum << endl;
    return read_pos_ranksum.get_raw_score();
}


double
snp_pos_info::
get_raw_baseQ()
{
//  cout << baseq_ranksum << endl;
    return baseq_ranksum.get_raw_score();
}


void
snp_pos_info::
print_known_counts(std::ostream& os,
                   const int min_qscore) const
{
    std::array<unsigned,N_BASE> base_count;
    get_known_counts(base_count,min_qscore);

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
    os.unsetf(std::ios::fixed);
}

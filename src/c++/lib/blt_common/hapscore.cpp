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
#include "blt_common/hapscore.hh"
#include "blt_util/math_util.hh"
#include "blt_util/qscore.hh"

#include <cmath>

#include <algorithm>
#include <iostream>
#include <vector>



hap_cand::
hap_cand(const bam_seq_base& read_seq,
         const uint8_t* init_qual,
         const int offset)  // the offset into read of the pileup base
    : _total_qual(0)
{
    const int read_len(read_seq.size());
    assert((offset>=0) && (offset<read_len));

    int start(offset-FLANK_SIZE);
    int end(offset+FLANK_SIZE+1);
    const int pre_seq( (start<0) ? -start : 0 );
    const int post_seq( (end>read_len) ? (end-read_len) : 0 );
    start=std::max(start,0);
    end=std::min(end,read_len);

    for (int i(0); i<pre_seq; ++i)
    {
        _bq[i] = 0;
    }

    for (int i(start); i<end; ++i)
    {
        _total_qual += init_qual[i];
        const char rs(read_seq.get_char(i));
        _bq[i-start+pre_seq] =
            ( (rs=='N') ?
              0 :
              (init_qual[i]<<QUAL_SHIFT | base_to_id(rs)));
    }

    for (int i(0); i<post_seq; ++i)
    {
        _bq[i+end-start+pre_seq] = 0;
    }
}


std::ostream&
operator<<(std::ostream& os, const hap_cand& hc)
{

    os << "hap_cand:\n";
    os << "total_qual: " << hc.total_qual() << "\n";
    os << "read: ";
    for (unsigned i(0); i<hap_cand::HAP_SIZE; ++i)
    {
        os << id_to_base(hc.base_id(i));
    }
    os << "\n";
    os << "qual: ";
    for (unsigned i(0); i<hap_cand::HAP_SIZE; ++i)
    {
        os << static_cast<char>(33+hc.qual(i));
    }
    os << "\n";
    return os;
}



struct hinfo
{
    hinfo() : total_qual(0) {}

    hinfo(const hap_cand& hc)
    {
        total_qual = hc.total_qual();
        for (unsigned i(0); i<hap_cand::HAP_SIZE; ++i)
        {
            hseq[i] = hc.base_id(i);
        }
    }

    bool
    operator<(const hinfo& rhs) const
    {
        return (total_qual > rhs.total_qual);
    }

    typedef boost::array<uint8_t,hap_cand::HAP_SIZE> hseq_t;
    hseq_t hseq;
    unsigned total_qual;
};


std::ostream&
operator<<(std::ostream& os, const hinfo& hi);



std::ostream&
operator<<(std::ostream& os, const hinfo& hi)
{

    os << "hinfo:\n";
    os << "total_qual: " << hi.total_qual << "\n";
    os << "read: ";
    for (unsigned i(0); i<hap_cand::HAP_SIZE; ++i)
    {
        os << id_to_base(hi.hseq[i]);
    }
    os << "\n";
    return os;
}



// tests for haplotype match
//
// if found (1) fills in any missing haplotype bases and (2) updates hinfo count
//
static
bool
is_hap_match(const hap_cand& hc,
             hinfo& hi)
{

    for (unsigned i(0); i<hap_cand::HAP_SIZE; ++i)
    {
        if (hi.hseq[i] != hc.base_id(i))
        {
            if ((hi.hseq[i] != BASE_ID::ANY) and
                (hc.base_id(i) != BASE_ID::ANY)) return false;
        }
    }

    // match!
    for (unsigned i(0); i<hap_cand::HAP_SIZE; ++i)
    {
        if ((hi.hseq[i] == hc.base_id(i)) or
            (hi.hseq[i] != BASE_ID::ANY)) continue;
        hi.hseq[i] = hc.base_id(i);
    }
    hi.total_qual += hc.total_qual();

    return true;
}



// alignment score is log(P(S|hap))
static
double
get_align_score(const hap_cand& hc,
                const hinfo& hi)
{

    //static const double lnany(std::log(0.25));
    static const double ln_one_third(-std::log(3.0));

    double score(0);
    for (unsigned i(0); i<hap_cand::HAP_SIZE; ++i)
    {
        if ((hi.hseq[i] == BASE_ID::ANY) ||
            (hc.base_id(i) == BASE_ID::ANY))
        {
            //score += lnany;
        }
        else if (hi.hseq[i] != hc.base_id(i))
        {
            score += qphred_to_ln_error_prob(hc.qual(i)) + ln_one_third - qphred_to_ln_comp_error_prob(hc.qual(i));
        }
    }

    return score;
}



double
get_hapscore(hap_set_t& hap_set)
{

    std::sort(hap_set.begin(),hap_set.end());

    std::vector<hinfo> haps;

    typedef hap_set_t::const_iterator hiter;
    {
        hiter i(hap_set.begin()), i_end(hap_set.end());
        for (; i!=i_end; ++i)
        {
            // 1: check if we match any types; add new type if not
            bool is_match(false);
            const unsigned hs(haps.size());
            for (unsigned j(0); j<hs; ++j)
            {
                if (is_hap_match(*i,haps[j]))
                {
                    is_match=true;
                    break;
                }
            }
            if (! is_match) haps.push_back(hinfo(*i));
        }
    }

    const bool is_2hap(haps.size()>1);

    // 2: get two most likely haplotypes:
    if (is_2hap)
    {
        std::partial_sort(haps.begin(),haps.begin()+2,haps.end());
    }

    // 3: calculate average read alignment score restricted to two best haplotypes:
    //static const double neginf(std::log(0));
    double ln2hapratio(0);
    {
        hiter i(hap_set.begin()), i_end(hap_set.end());
        for (; i!=i_end; ++i)
        {
            double als(get_align_score(*i,haps[0]));
            if (is_2hap)
            {
                als=std::max(als,get_align_score(*i,haps[1]));
            }
            //ln2hapratio = log_sum(ln2hapratio,als);
            ln2hapratio += als;
        }

        // target is log(avg(read_prob_ratio)), not avg(log(read_prob_ratio)):
        //
        //ln2hapratio -= std::log(static_cast<double>(hap_set.size()));
        ln2hapratio /= static_cast<double>(hap_set.size());
    }

    return ln2hapratio;
}

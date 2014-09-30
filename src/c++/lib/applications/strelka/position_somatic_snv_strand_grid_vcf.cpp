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

/// \author Chris Saunders
///

#include "position_somatic_snv_strand_grid_vcf.hh"
#include "strelka_vcf_locus_info.hh"
#include "somatic_call_shared.hh"
#if 0

#include "blt_common/snp_util.hh"
#include "blt_util/log.hh"
#include "blt_util/math_util.hh"
#include "blt_util/prob_util.hh"
#include "blt_util/seq_util.hh"

#include <cassert>
#include <cmath>
#include <cstdlib>

#include <map>

//#define SOMATIC_DEBUG

constexpr blt_float_t one_third(1./3.);
static const blt_float_t ln_one_third(std::log(one_third));
constexpr blt_float_t one_half(1./2.);
static const blt_float_t ln_one_half(std::log(one_half));
#endif
#include <fstream>
#include <iomanip>
#include <iostream>



static
void
write_vcf_sample_info(const blt_options& opt,
                      const extended_pos_data& tier1_epd,
                      const extended_pos_data& tier2_epd,
                      std::ostream& os)
{
    os << tier1_epd.n_calls
       << ':'
       << tier1_epd.n_unused_calls
       << ':'
       << tier1_epd.pi.n_spandel
       << ':'
       << tier1_epd.pi.n_submapped;

    std::array<unsigned,N_BASE> tier1_base_counts;
    std::array<unsigned,N_BASE> tier2_base_counts;
    tier1_epd.epd.good_pi.get_known_counts(tier1_base_counts,opt.used_allele_count_min_qscore);
    tier2_epd.epd.good_pi.get_known_counts(tier2_base_counts,opt.used_allele_count_min_qscore);
    for (unsigned b(0); b<N_BASE; ++b)
    {
        os << ':'
           << tier1_base_counts[b] << ','
           << tier2_base_counts[b];
    }
}



/// returns median, partially reorders elements in specified range
///
template <typename Iter>
typename std::iterator_traits<Iter>::value_type
median(
    Iter begin,
    Iter end)
{
    assert(begin != end);
    const auto size(std::distance(begin,end));
    std::nth_element(begin,begin+size/2, end);
    return *(begin+size/2);
}



static
double
safeFrac(const unsigned num, const unsigned denom)
{
    return ( (denom > 0) ? (num/static_cast<double>(denom)) : 0.);
}



void
write_vcf_somatic_snv_genotype_strand_grid(
    const strelka_options& opt,
    const strelka_deriv_options& dopt,
    const somatic_snv_genotype_grid& sgt,
    const bool is_write_nqss,
    const extended_pos_data& n1_epd,
    const extended_pos_data& t1_epd,
    const extended_pos_data& n2_epd,
    const extended_pos_data& t2_epd,
    std::ostream& os)
{
    typedef somatic_snv_genotype_grid::result_set result_set;

    const result_set& rs(sgt.rs);

    strelka_shared_modifiers smod;
    {
        // compute all site filters:
        const unsigned normalDP(n1_epd.n_calls);
        const unsigned tumorDP(t1_epd.n_calls);

        if (dopt.sfilter.is_max_depth())
        {
            if (normalDP > dopt.sfilter.max_depth)
            {
                smod.set_filter(STRELKA_VCF_FILTERS::HighDepth);
            }
        }

        {
            const unsigned normalFDP(n1_epd.n_unused_calls);
            const unsigned tumorFDP(t1_epd.n_unused_calls);

            const double normalFilt(safeFrac(normalFDP,normalDP));
            const double tumorFilt(safeFrac(tumorFDP,tumorDP));

            if ((normalFilt >=opt.sfilter.snv_max_filtered_basecall_frac) ||
                (tumorFilt >=opt.sfilter.snv_max_filtered_basecall_frac))
            {
                smod.set_filter(STRELKA_VCF_FILTERS::BCNoise);
            }
        }

        {
            const unsigned normalSDP(n1_epd.pi.n_spandel);
            const unsigned tumorSDP(t1_epd.pi.n_spandel);
            const unsigned normalSpanTot(normalDP + normalSDP);
            const unsigned tumorSpanTot(tumorDP + tumorSDP);

            const double normalSpanDelFrac(safeFrac(normalSDP, normalSpanTot));
            const double tumorSpanDelFrac(safeFrac(tumorSDP, tumorSpanTot));

            if ((normalSpanDelFrac > opt.sfilter.snv_max_spanning_deletion_frac) ||
                (tumorSpanDelFrac > opt.sfilter.snv_max_spanning_deletion_frac))
            {
                smod.set_filter(STRELKA_VCF_FILTERS::SpanDel);
            }
        }

        if ((rs.ntype != NTYPE::REF) || (rs.snv_from_ntype_qphred < opt.sfilter.snv_min_qss_ref))
        {
            smod.set_filter(STRELKA_VCF_FILTERS::QSS_ref);
        }
    }

    //REF:
    os << '\t' << n1_epd.pi.ref_base
       //ALT:
       << "\t";
    DDIGT_SGRID::write_alt_alleles(static_cast<DDIGT_SGRID::index_t>(rs.max_gt),
                                   sgt.ref_gt,os);
    //QUAL:
    os << "\t.";

    //FILTER:
    os << "\t";
    smod.write_filters(os);

    //INFO:
    os << '\t'
       << "SOMATIC"
       << ";QSS=" << rs.snv_qphred;

    if (is_write_nqss)
    {
        os << ";NQSS=" << rs.nonsomatic_qphred;
    }

    os << ";TQSS=" << (sgt.snv_tier+1)
       << ";NT=" << NTYPE::label(rs.ntype)
       << ";QSS_NT=" << rs.snv_from_ntype_qphred
       << ";TQSS_NT=" << (sgt.snv_from_ntype_tier+1)
       << ";SGT=";
    DDIGT_SGRID::write_state(static_cast<DDIGT_SGRID::index_t>(rs.max_gt),
                             sgt.ref_gt,os);

    {
        std::ofstream tmp_os;
        tmp_os.copyfmt(os);
        os << std::fixed << std::setprecision(2);

        {
            // m_mapq includes all calls, even from reads below the mapq threshold:
            const unsigned n_mapq(n1_epd.pi.n_mapq+t1_epd.pi.n_mapq);
            os << ";DP=" << n_mapq;
        }

        {
            const unsigned n_mapq(n1_epd.pi.n_mapq+t1_epd.pi.n_mapq);
            const double cumm_mapq2(n1_epd.pi.cumm_mapq + t1_epd.pi.cumm_mapq);

            os << ";MQ=" << std::sqrt(cumm_mapq2/n_mapq);
        }

        {
            const unsigned n_mapq0(n1_epd.pi.n_mapq0+t1_epd.pi.n_mapq0);

            os << ";MQ0=" << n_mapq0;
        }

        {
            // don't worry about efficiency for median calc right now:
            bool isAltpos(false);
            bool isAltmap(false);
            uint16_t altpos=0;
            uint16_t altmap=0;
            if (! t1_epd.pi.altReadPos.empty())
            {
                const auto& apos(t1_epd.pi.altReadPos);
                std::vector<uint16_t> readpos;
                for (const auto& r : apos)
                {
                    readpos.push_back(r.readPos);
                }
                std::vector<uint16_t> readposcomp;
                for (const auto& r : apos)
                {
                    readposcomp.push_back(r.readPosLength-r.readPos);
                }

                const auto pmedian(median(readpos.begin(),readpos.end()));
                const auto lmedian(median(readposcomp.begin(),readposcomp.end()));

                isAltpos=true;
                altpos=std::min(pmedian,lmedian);

                if (readpos.size() >= 3)
                {
                    for (auto& p : readpos)
                    {
                        p = std::abs(p-pmedian);
                    }

                    isAltmap=true;
                    altmap=median(readpos.begin(),readpos.end());
                }
            }

            os << ";ALTPOS=";
            if (isAltpos) os << altpos;
            else          os << '.';

            os << ";ALTMAP=";
            if (isAltmap) os << altmap;
            else          os << '.';
        }

        {
            const double ReadPosRankSum = t1_epd.pi.read_pos_ranksum.get_u_stat();
            os << ";ReadPosRankSum=" << ReadPosRankSum;
        }

        os.copyfmt(tmp_os);
    }

    //FORMAT:
    os << '\t'
       << "DP:FDP:SDP:SUBDP:AU:CU:GU:TU";

    // normal sample info:
    os << "\t";
    write_vcf_sample_info(opt,n1_epd,n2_epd,os);

    // tumor sample info:
    os << "\t";
    write_vcf_sample_info(opt,t1_epd,t2_epd,os);
}

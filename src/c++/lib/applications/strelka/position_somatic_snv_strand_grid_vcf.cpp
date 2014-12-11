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
/// \author Morten Kallberg
///

#include "position_somatic_snv_strand_grid_vcf.hh"
#include "strelka_vcf_locus_info.hh"
#include "somatic_call_shared.hh"
#include "blt_util/io_util.hh"
#include "blt_util/math_util.hh"

#include <iomanip>
#include <iostream>


template <typename D>
static
double
safeFrac(const unsigned num, const D denom)
{
    return ( (denom > 0) ? (num/static_cast<double>(denom)) : 0.);
}



// similar to 'write_vcf_sample_info' below, redundancy needed to get order of output right
// TODO consolidate this dual calculation step
static
void
set_VQSR_sample_info(
    const blt_options& opt,
    const strelka_deriv_options& dopt,
    const extended_pos_data& tier1_epd,
    const extended_pos_data& tier2_epd,
    strelka_shared_modifiers& smod,
    const char ref_base,
    bool isNormalSample)
{
    using namespace STRELKA_VQSR_FEATURES;

    // add in VQSR features
    // {2,N_FDP_RATE},{3,T_FDP_RATE},{4,N_SDP_RATE},
    // {5,T_SDP_RATE},{6,N_DP_RATE},{7,TIER1_ALLELE_RATE}
    const double FDP_ratio(safeFrac(tier1_epd.n_unused_calls, tier1_epd.n_calls));
    const double SDP_ratio(safeFrac(tier1_epd.pi.n_spandel, tier1_epd.n_calls+tier1_epd.pi.n_spandel));

    smod.set_feature((isNormalSample ? N_FDP_RATE : T_FDP_RATE), FDP_ratio);
    smod.set_feature((isNormalSample ? N_SDP_RATE : T_SDP_RATE), SDP_ratio);

    if (isNormalSample)      // offset of 1 is tumor case, we only calculate the depth rate for the normal
    {
        smod.set_feature(N_DP_RATE, safeFrac(tier1_epd.n_calls,dopt.sfilter.max_depth));
    }

    if (!isNormalSample)      //report tier1_allele count for tumor case
    {
        std::array<unsigned,N_BASE> tier1_base_counts;
        tier1_epd.epd.good_pi.get_known_counts(tier1_base_counts,opt.used_allele_count_min_qscore);

        const unsigned ref_index(base_to_id(ref_base));
        unsigned ref=0;
        unsigned alt=0;
        for (unsigned b(0); b<N_BASE; ++b)
        {
            if (b==ref_index)
            {
                ref += tier1_base_counts[b];
            }
            else
            {
                alt += tier1_base_counts[b];
            }
        }

        const double allele_freq(safeFrac(alt,ref+alt));
        smod.set_feature(TIER1_ALLELE_RATE,allele_freq);
    }
}



// Prepare feature vector in case we are using VQSR, the individual values will be set
// feature_type ft;

// Feature dictionary constructed position-> feature similar to the columns in the training set (correct ordering is key)
// QSS_NT, N_FDP_RATE,T_FDP_RATE,N_SDP_RATE,T_SDP_RATE,N_DP_RATE,TIER1_ALLELE_RATE,MQ_SCORE
// MQ_ZERO_RATE,SNVSB,ReadPosRankSum,ALTMAP,ALTPOS,PNOISE,PNOISE2
//    Definition of features that need calculated
//    # Compute the FDP, SDP, and DP rates.
//    n_FDP_ratio     = n_FDP/n_DP if n_DP != 0 else 0
//    t_FDP_ratio     = t_FDP/t_DP if t_DP != 0 else 0
//
//    n_SDP_ratio     = n_SDP/(n_DP + n_SDP)  if (n_DP + n_SDP) != 0 else 0
//    t_SDP_ratio     = t_SDP/(t_DP + t_SDP ) if (t_DP + t_SDP) != 0 else 0
//    n_DP_ratio      = n_DP/float(self.chr_depth)
//    t_DP_ratio      = t_DP/float(self.chr_depth)
//    mq_zero_ratio   = MQ_ZERO/float(self.chr_depth)
//
//    Allele count logic
//    if t_allele_alt_counts[0] + t_allele_ref_counts[0] == 0:
//        t_tier1_allele_rate = 0
//    else:
//        t_tier1_allele_rate = t_allele_alt_counts[0] / float(t_allele_alt_counts[0] + t_allele_ref_counts[0])
//    ft = {{1,rs.snv_from_ntype_qphred},{2,N_FDP_RATE},{3,T_FDP_RATE},{4,N_SDP_RATE},
//          {5,T_SDP_RATE},{6,N_DP_RATE},{7,TIER1_ALLELE_RATE},{8,MQ},
//          {9,n_mapq0},{10,rs.strandBias},{11,ReadPosRankSum},{12,altmap},{13,altpos},{14,pnoise},{15,pnoise2}};
static
void
calc_VQSR_features(
    const blt_options& opt,
    const strelka_deriv_options& dopt,
    const somatic_snv_genotype_grid& sgt,
    strelka_shared_modifiers& smod,
    const extended_pos_data& n1_epd,
    const extended_pos_data& t1_epd,
    const extended_pos_data& n2_epd,
    const extended_pos_data& t2_epd,
    const result_set& rs)
{
    // don't worry about efficiency for median calc right now:
    //bool isAltpos(false);
    //bool isAltmap(false);
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

        // isAltpos=true;
        altpos=std::min(pmedian,lmedian);

        if (readpos.size() >= 3)
        {
            for (auto& p : readpos)
            {
                p = std::abs(p-pmedian);
            }

            //isAltmap=true;
            altmap=median(readpos.begin(),readpos.end());
        }
    }

    //QSS_NT
    smod.set_feature(STRELKA_VQSR_FEATURES::QSS_NT,rs.snv_from_ntype_qphred);

    set_VQSR_sample_info(opt,dopt,n1_epd,n2_epd,smod,n1_epd.pi.ref_base,true);
    set_VQSR_sample_info(opt,dopt,t1_epd,t2_epd,smod,n1_epd.pi.ref_base,false);

    //MQ
    const unsigned n_mapq(n1_epd.pi.n_mapq+t1_epd.pi.n_mapq);
    const double cumm_mapq2(n1_epd.pi.cumm_mapq + t1_epd.pi.cumm_mapq);
    smod.set_feature(STRELKA_VQSR_FEATURES::MQ,std::sqrt(cumm_mapq2/n_mapq));

    //n_mapq0
    const unsigned n_mapq0(n1_epd.pi.n_mapq0+t1_epd.pi.n_mapq0);
    smod.set_feature(STRELKA_VQSR_FEATURES::n_mapq0, safeFrac(n_mapq0,n_mapq0+n_mapq));

    //ReadPosRankSum
    const double ReadPosRankSum = t1_epd.pi.read_pos_ranksum.get_u_stat();
    smod.set_feature(STRELKA_VQSR_FEATURES::ReadPosRankSum,ReadPosRankSum);

    //StrandBias
    smod.set_feature(STRELKA_VQSR_FEATURES::strandBias,rs.strandBias);

    /// TODO better handling of default values for in cases where alpos or altmap are not defined (0 is not a good default)
    ///
    //Altpos
    smod.set_feature(STRELKA_VQSR_FEATURES::altpos,altpos);
    smod.set_feature(STRELKA_VQSR_FEATURES::altmap,altmap);

    //Pnoise
    double pnoise(0);
    if (sgt.sn.total > 1 && sgt.sn.noise > 1)
    {
        pnoise = sgt.sn.nfrac();
    }
    smod.set_feature(STRELKA_VQSR_FEATURES::pnoise,pnoise);

    //Pnoise2
    double pnoise2(0);
    if (sgt.sn.total > 1 && sgt.sn.noise2 > 1)
    {
        pnoise2 = sgt.sn.n2frac();
    }
    smod.set_feature(STRELKA_VQSR_FEATURES::pnoise2,pnoise2);
}



static
void
write_vcf_sample_info(
    const blt_options& opt,
    const strelka_deriv_options& dopt,
    const extended_pos_data& tier1_epd,
    const extended_pos_data& tier2_epd,
    std::ostream& os)
{
    //DP:FDP:SDP:SUBDP:AU:CU:GU:TU
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

    {
        // Make sure the VQSR feature vector is populated
        // this is done even if not running with VQSR as some intermediate
        // calculations are still needed for VCF reporting
        calc_VQSR_features(opt,dopt,sgt,smod,n1_epd,t1_epd,n2_epd,t2_epd,rs);

        // case we are doing VQSR, clear filters and apply single LowQscore filter
        if (scoring_models::Instance()->calibration_init)  // write out somatic VQSR metrics
        {
            smod.Qscore = scoring_models::Instance()->score_instance(smod.get_features());
            smod.filters.reset();

            // Temp hack to handle sample with large LOH, if REF is already het, set low score and filter by default
            if (rs.ntype != NTYPE::REF) smod.Qscore=0;

            if (smod.Qscore < opt.sfilter.minimumQscore)
                smod.set_filter(STRELKA_VCF_FILTERS::LowQscore);
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
       << ";TQSS_NT=" << (sgt.snv_from_ntype_tier+1);
    os << ";SGT=";

    DDIGT_SGRID::write_state(static_cast<DDIGT_SGRID::index_t>(rs.max_gt),
                             sgt.ref_gt,os);

    {
        const StreamScoper ss(os);
        os << std::fixed << std::setprecision(2);

        // m_mapq includes all calls, even from reads below the mapq threshold:
        const unsigned n_mapq(n1_epd.pi.n_mapq+t1_epd.pi.n_mapq);
        os << ";DP=" << n_mapq;
        os << ";MQ=" << smod.get_feature(STRELKA_VQSR_FEATURES::MQ);

        {
            const unsigned n_mapq0(n1_epd.pi.n_mapq0+t1_epd.pi.n_mapq0);
            os << ";MQ0=" << n_mapq0;
        }

        os << ";ALTPOS=";
        if (smod.test_feature(STRELKA_VQSR_FEATURES::altpos))
            os << (int)smod.get_feature(STRELKA_VQSR_FEATURES::altpos);
        else
            os << '.';

        os << ";ALTMAP=";
        if (smod.test_feature(STRELKA_VQSR_FEATURES::altmap))
            os << (int)smod.get_feature(STRELKA_VQSR_FEATURES::altmap);
        else
            os << '.';

        os << ";ReadPosRankSum=" << smod.get_feature(STRELKA_VQSR_FEATURES::ReadPosRankSum);
        os << ";SNVSB=" << smod.get_feature(STRELKA_VQSR_FEATURES::strandBias);
        os << ";PNOISE=" << smod.get_feature(STRELKA_VQSR_FEATURES::pnoise);
        os << ";PNOISE2=" << smod.get_feature(STRELKA_VQSR_FEATURES::pnoise2);

        if (scoring_models::Instance()->calibration_init)
        {
            os << ";VQSR=" << smod.Qscore;
        }
    }

    //FORMAT:
    os << '\t'
       << "DP:FDP:SDP:SUBDP:AU:CU:GU:TU";

    // normal sample info:
    os << "\t";
    write_vcf_sample_info(opt,dopt,n1_epd,n2_epd,os);

    // tumor sample info:
    os << "\t";
    write_vcf_sample_info(opt,dopt,t1_epd,t2_epd,os);
}

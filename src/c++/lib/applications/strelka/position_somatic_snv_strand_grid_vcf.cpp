// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2016 Illumina, Inc.
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

/// \author Chris Saunders
/// \author Morten Kallberg
///

#include "position_somatic_snv_strand_grid_vcf.hh"
#include "somaticAlleleUtil.hh"
#include "strelka_vcf_locus_info.hh"
#include "somatic_call_shared.hh"
#include "blt_util/io_util.hh"
#include "blt_util/math_util.hh"


#include <iomanip>
#include <iostream>



double
calculateLogOddsRatio(
    const CleanedPileup& n1_cpi,
    const CleanedPileup& t1_cpi,
    const blt_options& opt)
{
    /**
     * Calculate LOR feature (log odds ratio for  T_REF T_ALT
     *                                            N_REF N_ALT)
     */
    std::array<unsigned,N_BASE> n_tier1_base_counts;
    std::array<unsigned,N_BASE> t_tier1_base_counts;
    n1_cpi.cleanedPileup().get_known_counts(n_tier1_base_counts,opt.used_allele_count_min_qscore);
    t1_cpi.cleanedPileup().get_known_counts(t_tier1_base_counts,opt.used_allele_count_min_qscore);
    const unsigned ref_index(base_to_id(n1_cpi.cleanedPileup().get_ref_base()));
    unsigned n_ref=0;
    unsigned n_alt=0;
    unsigned t_ref=0;
    unsigned t_alt=0;
    for (unsigned b(0); b<N_BASE; ++b)
    {
        if (b==ref_index)
            n_ref += n_tier1_base_counts[b];
        else
            n_alt += n_tier1_base_counts[b];
    }
    for (unsigned b(0); b<N_BASE; ++b)
    {
        if (b==ref_index)
            t_ref += t_tier1_base_counts[b];
        else
            t_alt += t_tier1_base_counts[b];
    }
    double LOR = log10((t_ref+0.5)*(n_alt+0.5) / (t_alt+0.5) / (n_ref+0.5));

    return LOR;
}



/// set sample specific empirical scoring features
///
// similar to 'write_vcf_sample_info' below, redundancy needed to get order of output right
// TODO consolidate this dual calculation step
static
void
get_single_sample_scoring_features(
    const blt_options& opt,
    const strelka_deriv_options& dopt,
    const CleanedPileup& tier1_cpi,
    const CleanedPileup& /*tier2_cpi*/,
    const double normChromDepth,
    const bool isNormalSample,
    strelka_shared_modifiers_snv& smod)
{
    if (opt.isReportEVSFeatures)
    {
        // {2,N_FDP_RATE},{3,T_FDP_RATE},{4,N_SDP_RATE},
        // {5,T_SDP_RATE},{6,N_DP_RATE},{7,TIER1_ALLELE_RATE}
        const double FDP_ratio(safeFrac(tier1_cpi.n_unused_calls(), tier1_cpi.n_calls()));
        const double SDP_ratio(safeFrac(tier1_cpi.rawPileup().spanningDeletionReadCount, tier1_cpi.n_calls()+tier1_cpi.rawPileup().spanningDeletionReadCount));

        smod.dfeatures.set((isNormalSample ? SOMATIC_SNV_SCORING_DEVELOPMENT_FEATURES::N_FDP_RATE : SOMATIC_SNV_SCORING_DEVELOPMENT_FEATURES::T_FDP_RATE), FDP_ratio);
        smod.dfeatures.set((isNormalSample ? SOMATIC_SNV_SCORING_DEVELOPMENT_FEATURES::N_SDP_RATE : SOMATIC_SNV_SCORING_DEVELOPMENT_FEATURES::T_SDP_RATE), SDP_ratio);
    }

    // compute N_DP_RATE
    const bool isUniformDepthExpected(dopt.sfilter.is_max_depth());
    if (isNormalSample)      // offset of 1 is tumor case, we only calculate the depth rate for the normal
    {
        double normalDepthRate(1.);
        if (isUniformDepthExpected)
        {
            normalDepthRate = safeFrac(tier1_cpi.n_calls(), normChromDepth);
        }

        smod.features.set(SOMATIC_SNV_SCORING_FEATURES::N_DP_RATE,normalDepthRate);
    }


    if (!isNormalSample)      //report tier1_allele count for tumor case
    {
        std::array<unsigned,N_BASE> tier1_base_counts;
        tier1_cpi.cleanedPileup().get_known_counts(tier1_base_counts,opt.used_allele_count_min_qscore);

        const unsigned ref_index(base_to_id(tier1_cpi.cleanedPileup().get_ref_base()));
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

        // cap the allele rate at 0.5
        smod.features.set(SOMATIC_SNV_SCORING_FEATURES::TIER1_ALT_RATE,std::min(0.5,allele_freq));
    }
}



// Prepare feature vector in case we are using empirical scoring, the individual values will be set
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
get_scoring_features(
    const blt_options& opt,
    const strelka_deriv_options& dopt,
    const somatic_snv_genotype_grid& /*sgt*/,
    const CleanedPileup& n1_cpi,
    const CleanedPileup& t1_cpi,
    const CleanedPileup& n2_cpi,
    const CleanedPileup& t2_cpi,
    const double normChromDepth,
    const snv_result_set& rs,
    strelka_shared_modifiers_snv& smod)
{
    uint16_t altpos=0;
    uint16_t altmap=0;
    if (! t1_cpi.rawPileup().altReadPos.empty())
    {
        const auto& apos(t1_cpi.rawPileup().altReadPos);
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

        altpos=std::min(pmedian,lmedian);

        if (readpos.size() >= 3)
        {
            for (auto& p : readpos)
            {
                p = std::abs(p-pmedian);
            }

            altmap=median(readpos.begin(),readpos.end());
        }
    }

    //QSS_NT
    smod.features.set(SOMATIC_SNV_SCORING_FEATURES::QSS_NT,rs.from_ntype_qphred);

    static const bool isNormalSample(true);
    get_single_sample_scoring_features(opt,dopt,n1_cpi,n2_cpi, normChromDepth, isNormalSample,smod);
    get_single_sample_scoring_features(opt,dopt,t1_cpi,t2_cpi, normChromDepth, (!isNormalSample),smod);

    //MQ
    MapqTracker mapqTracker(n1_cpi.rawPileup().mapqTracker);
    mapqTracker.merge(t1_cpi.rawPileup().mapqTracker);
    smod.features.set(SOMATIC_SNV_SCORING_FEATURES::MQ, mapqTracker.getRMS());

    //n_mapq0
    const unsigned n_mapq(mapqTracker.count);
    const unsigned n_mapq0(mapqTracker.zeroCount);
    smod.features.set(SOMATIC_SNV_SCORING_FEATURES::MQ0_FRAC, safeFrac(n_mapq0,n_mapq));

    //ReadPosRankSum
    const double ReadPosRankSum = t1_cpi.rawPileup().read_pos_ranksum.get_z_stat();
    smod.features.set(SOMATIC_SNV_SCORING_FEATURES::ReadPosRankSum,ReadPosRankSum);

    //StrandBias
    smod.features.set(SOMATIC_SNV_SCORING_FEATURES::strandBias,rs.strandBias);

    smod.features.set(SOMATIC_SNV_SCORING_FEATURES::LOR, calculateLogOddsRatio(n1_cpi,t1_cpi,opt));

    // features not used in the current EVS model but feature candidates/exploratory for new EVS models
    if (opt.isReportEVSFeatures)
    {
        /// TODO better handling of default values for in cases where altpos or altmap are not defined (0 is not a good default)
        smod.dfeatures.set(SOMATIC_SNV_SCORING_DEVELOPMENT_FEATURES::altpos,altpos);
        smod.dfeatures.set(SOMATIC_SNV_SCORING_DEVELOPMENT_FEATURES::altmap,altmap);
    }
}


static
void
write_vcf_sample_info(
    const blt_options& opt,
    const strelka_deriv_options& /*dopt*/,
    const CleanedPileup& tier1_cpi,
    const CleanedPileup& tier2_cpi,
    std::ostream& os)
{
    //DP:FDP:SDP:SUBDP:AU:CU:GU:TU
    os << tier1_cpi.n_calls()
       << ':'
       << tier1_cpi.n_unused_calls()
       << ':'
       << tier1_cpi.rawPileup().spanningDeletionReadCount
       << ':'
       << tier1_cpi.rawPileup().n_submapped;


    std::array<unsigned,N_BASE> tier1_base_counts;
    std::array<unsigned,N_BASE> tier2_base_counts;
    tier1_cpi.cleanedPileup().get_known_counts(tier1_base_counts,opt.used_allele_count_min_qscore);
    tier2_cpi.cleanedPileup().get_known_counts(tier2_base_counts,opt.used_allele_count_min_qscore);

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
    const CleanedPileup& n1_epd,
    const CleanedPileup& t1_epd,
    const CleanedPileup& n2_epd,
    const CleanedPileup& t2_epd,
    const double normChromDepth,
    const double maxChromDepth,
    std::ostream& os)
{
    const snv_result_set& rs(sgt.rs);

    strelka_shared_modifiers_snv smod;

    const bool isEVS(dopt.somaticSnvScoringModel);

    if (! isEVS)
    {
        // compute all site filters:
        const unsigned normalDP(n1_epd.n_calls());
        const unsigned tumorDP(t1_epd.n_calls());

        if (dopt.sfilter.is_max_depth())
        {
            if (normalDP > maxChromDepth)
            {
                smod.filters.set(SOMATIC_VARIANT_VCF_FILTERS::HighDepth);
            }
        }

        {
            const unsigned normalFDP(n1_epd.n_unused_calls());
            const unsigned tumorFDP(t1_epd.n_unused_calls());

            const double normalFilt(safeFrac(normalFDP,normalDP));
            const double tumorFilt(safeFrac(tumorFDP,tumorDP));

            if ((normalFilt >=opt.sfilter.snv_max_filtered_basecall_frac) ||
                (tumorFilt >=opt.sfilter.snv_max_filtered_basecall_frac))
            {
                smod.filters.set(SOMATIC_VARIANT_VCF_FILTERS::BCNoise);
            }
        }

        {
            const unsigned normalSDP(n1_epd.rawPileup().spanningDeletionReadCount);
            const unsigned tumorSDP(t1_epd.rawPileup().spanningDeletionReadCount);
            const unsigned normalSpanTot(normalDP + normalSDP);
            const unsigned tumorSpanTot(tumorDP + tumorSDP);

            const double normalSpanDelFrac(safeFrac(normalSDP, normalSpanTot));
            const double tumorSpanDelFrac(safeFrac(tumorSDP, tumorSpanTot));

            if ((normalSpanDelFrac > opt.sfilter.snv_max_spanning_deletion_frac) ||
                (tumorSpanDelFrac > opt.sfilter.snv_max_spanning_deletion_frac))
            {
                smod.filters.set(SOMATIC_VARIANT_VCF_FILTERS::SpanDel);
            }
        }

        if ((rs.ntype != NTYPE::REF) || (rs.from_ntype_qphred < opt.sfilter.snv_min_qss_ref))
        {
            smod.filters.set(SOMATIC_VARIANT_VCF_FILTERS::QSS_ref);
        }
    }

    {
        // Make sure the empirical scoring feature vector is populated
        // this is done even if not running with EVS as some intermediate
        // calculations are still needed for VCF reporting
        get_scoring_features(opt,dopt,sgt,n1_epd,t1_epd,n2_epd,t2_epd, normChromDepth, rs,smod);

        // if we are using empirical scoring, clear filters and apply single LowEVS filter
        if (isEVS)
        {
            assert(dopt.somaticSnvScoringModel);
            const VariantScoringModelServer& varModel(*dopt.somaticSnvScoringModel);
            updateAlleleEVSScore(varModel, rs, smod);

            smod.filters.clear();
            if (smod.EVS < varModel.scoreFilterThreshold())
            {
                smod.filters.set(SOMATIC_VARIANT_VCF_FILTERS::LowEVSsnv);
            }
        }
    }

    char ref_base = n1_epd.rawPileup().get_ref_base();
    //REF:
    os << '\t' << ref_base
       //ALT:
       << "\t";

    DDIGT::write_alt_alleles(id_to_base(rs.normal_alt_id),
                             id_to_base(rs.tumor_alt_id),
                             ref_base,
                             os);

    //QUAL:
    os << "\t.";

    //FILTER:
    os << "\t";
    smod.filters.write(os);

    //INFO:
    os << '\t'
       << "SOMATIC"
       << ";QSS=" << rs.qphred;

    if (is_write_nqss)
    {
        os << ";NQSS=" << rs.nonsomatic_qphred;
    }

    os << ";TQSS=" << (sgt.snv_tier+1)
       << ";NT=" << NTYPE::label(rs.ntype)
       << ";QSS_NT=" << rs.from_ntype_qphred
       << ";TQSS_NT=" << (sgt.snv_from_ntype_tier+1);
    os << ";SGT=";

    DDIGT::write_snv_state(static_cast<DDIGT::index_t>(rs.max_gt),
                           ref_base,
                           id_to_base(rs.normal_alt_id),
                           id_to_base(rs.tumor_alt_id),
                           os);

    {
        const StreamScoper ss(os);
        os << std::fixed << std::setprecision(2);

        // m_mapq includes all calls, even from reads below the mapq threshold:
        MapqTracker mapqTracker(n1_epd.rawPileup().mapqTracker);
        mapqTracker.merge(t1_epd.rawPileup().mapqTracker);
        os << ";DP=" << mapqTracker.count;
        os << ";MQ=" << smod.features.get(SOMATIC_SNV_SCORING_FEATURES::MQ);
        os << ";MQ0=" << mapqTracker.zeroCount;

//        os << ";ALTPOS=";
//        if (smod.features.test(SOMATIC_SNV_SCORING_FEATURES::altpos))
//            os << (int)smod.features.get(SOMATIC_SNV_SCORING_FEATURES::altpos);
//        else
//            os << '.';

//        os << ";ALTMAP=";
//        if (smod.features.test(SOMATIC_SNV_SCORING_FEATURES::altmap))
//            os << (int)smod.features.get(SOMATIC_SNV_SCORING_FEATURES::altmap);
//        else
//            os << '.';

        os << ";ReadPosRankSum=" << smod.features.get(SOMATIC_SNV_SCORING_FEATURES::ReadPosRankSum);
        os << ";SNVSB=" << smod.features.get(SOMATIC_SNV_SCORING_FEATURES::strandBias);

        if (smod.isEVS)
        {
            os << ";EVS=" << smod.EVS;
        }
    }

    if (opt.isReportEVSFeatures)
    {
        const StreamScoper ss(os);
        os << std::setprecision(5);
        os << ";EVSF=";
        smod.features.writeValues(os);
        os << ",";
        smod.dfeatures.writeValues(os);
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

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


/// Calculate LOR feature for SNVs (log odds ratio for  T_REF T_ALT
///                               		                N_REF N_ALT)
///
static
double
calculateLogOddsRatio(
    const CleanedPileup& n1_cpi,
    const CleanedPileup& t1_cpi)
{
    std::array<unsigned,N_BASE> n_tier1_base_counts;
    std::array<unsigned,N_BASE> t_tier1_base_counts;
    n1_cpi.cleanedPileup().getBasecallCounts(n_tier1_base_counts);
    t1_cpi.cleanedPileup().getBasecallCounts(t_tier1_base_counts);
    const unsigned ref_index(base_to_id(n1_cpi.cleanedPileup().get_ref_base()));

    static const double pseudoCount(0.5);
    double normalRefCount(pseudoCount);
    double normalAltCount(pseudoCount);
    double tumorRefCount(pseudoCount);
    double tumorAltCount(pseudoCount);
    for (unsigned b(0); b<N_BASE; ++b)
    {
        if (b==ref_index)
        {
            normalRefCount += n_tier1_base_counts[b];
            tumorRefCount += t_tier1_base_counts[b];
        }
        else
        {
            normalAltCount += n_tier1_base_counts[b];
            tumorAltCount += t_tier1_base_counts[b];
        }
    }

    return std::log((tumorRefCount*normalAltCount) / (tumorAltCount*normalRefCount));
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
    const double filteredDepthFraction(safeFrac(tier1_cpi.unusedBasecallCount(), tier1_cpi.totalBasecallCount()));
    smod.features.set((isNormalSample ? SOMATIC_SNV_SCORING_FEATURES::NormalSampleFilteredDepthFraction : SOMATIC_SNV_SCORING_FEATURES::TumorSampleFilteredDepthFraction), filteredDepthFraction);
    if (opt.isReportEVSFeatures)
    {
        const double spanningDeletionFraction(safeFrac(tier1_cpi.rawPileup().spanningDeletionReadCount,
                                                       tier1_cpi.totalBasecallCount()+tier1_cpi.rawPileup().spanningDeletionReadCount));
        smod.dfeatures.set((isNormalSample ? SOMATIC_SNV_SCORING_DEVELOPMENT_FEATURES::NormalSampleSpanningDeletionFraction : SOMATIC_SNV_SCORING_DEVELOPMENT_FEATURES::TumorSampleSpanningDeletionFraction), spanningDeletionFraction);
    }

    // compute NormalSampleRelativeTotalLocusDepth
    const bool isUniformDepthExpected(dopt.sfilter.is_max_depth());
    if (isNormalSample)      // offset of 1 is tumor case, we only calculate the depth rate for the normal
    {
        double normalDepthRate(1.);
        if (isUniformDepthExpected)
        {
            /// TODO in theory we would like the numerator to include all reads that are used to compute normChromDepth -- I'm not sure this is the case now:
            normalDepthRate = safeFrac(tier1_cpi.totalBasecallCount(), normChromDepth);
        }

        smod.features.set(SOMATIC_SNV_SCORING_FEATURES::NormalSampleRelativeTotalLocusDepth,normalDepthRate);
    }


    if (!isNormalSample)      //report tier1_allele count for tumor case
    {
        std::array<unsigned,N_BASE> tier1_base_counts;
        tier1_cpi.cleanedPileup().getBasecallCounts(tier1_base_counts);

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

        // cap the allele rate at 0.5 to help prevent the EVS model from overtraining against LOH regions with
        // higher allele frequencies
        smod.features.set(SOMATIC_SNV_SCORING_FEATURES::TumorSampleAltAlleleFraction,std::min(0.5,allele_freq));
    }
}



/// Compute variant features needed for the scoring model and/or direct VCF output
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
    uint16_t medianReadPos=0;
    uint16_t medianReadPosVar=0;
    if (! t1_cpi.rawPileup().altAlleleReadPositionInfo.empty())
    {
        const auto& apos(t1_cpi.rawPileup().altAlleleReadPositionInfo);
        std::vector<uint16_t> readpos;
        for (const auto& r : apos)
        {
            readpos.push_back(r.readPos);
        }
        std::vector<uint16_t> readposcomp;
        for (const auto& r : apos)
        {
            readposcomp.push_back(r.readLength-r.readPos);
        }

        const auto pmedian(median(readpos.begin(),readpos.end()));
        const auto lmedian(median(readposcomp.begin(),readposcomp.end()));

        medianReadPos=std::min(pmedian,lmedian);

        if (readpos.size() >= 3)
        {
            for (auto& p : readpos)
            {
                p = std::abs(p-pmedian);
            }

            medianReadPosVar=median(readpos.begin(),readpos.end());
        }
    }

    {
        const int from_ref_qphred((rs.ntype == NTYPE::REF) ? rs.from_ntype_qphred : 0 );
        smod.features.set(SOMATIC_SNV_SCORING_FEATURES::SomaticSNVQualityAndHomRefGermlineGenotype, from_ref_qphred);
    }

    static const bool isNormalSample(true);
    get_single_sample_scoring_features(opt,dopt,n1_cpi,n2_cpi, normChromDepth, isNormalSample,smod);
    get_single_sample_scoring_features(opt,dopt,t1_cpi,t2_cpi, normChromDepth, (!isNormalSample),smod);

    MapqTracker mapqTracker(n1_cpi.rawPileup().mapqTracker);
    mapqTracker.merge(t1_cpi.rawPileup().mapqTracker);
    smod.features.set(SOMATIC_SNV_SCORING_FEATURES::RMSMappingQuality, mapqTracker.getRMS());

    const unsigned mappingQualityObservations(mapqTracker.count);
    const unsigned zeroMappingQualityObservations(mapqTracker.zeroCount);
    smod.features.set(SOMATIC_SNV_SCORING_FEATURES::ZeroMappingQualityFraction, safeFrac(zeroMappingQualityObservations,mappingQualityObservations));

    const double tumorSampleReadPosRankSum = t1_cpi.rawPileup().readPositionRankSum.get_z_stat();
    smod.features.set(SOMATIC_SNV_SCORING_FEATURES::TumorSampleReadPosRankSum, tumorSampleReadPosRankSum);

    smod.features.set(SOMATIC_SNV_SCORING_FEATURES::TumorSampleStrandBias,rs.strandBias);

    smod.features.set(SOMATIC_SNV_SCORING_FEATURES::AlleleCountLogOddsRatio, calculateLogOddsRatio(n1_cpi, t1_cpi));

    // features not used in the current EVS model but feature candidates/exploratory for new EVS models
    if (opt.isReportEVSFeatures)
    {
        smod.dfeatures.set(SOMATIC_SNV_SCORING_DEVELOPMENT_FEATURES::TumorSampleAltAlleleMedianReadPos,medianReadPos);
        smod.dfeatures.set(SOMATIC_SNV_SCORING_DEVELOPMENT_FEATURES::TumorSampleAltAlleleMedianReadPosVariation,medianReadPosVar);
    }
}



/// Write all values for the VCF's tumor or normal SAMPLE field.
static
void
write_vcf_sample_info(
    const strelka_deriv_options& /*dopt*/,
    const CleanedPileup& tier1_cpi,
    const CleanedPileup& tier2_cpi,
    std::ostream& os)
{
    //DP:FDP:SDP:SUBDP:AU:CU:GU:TU
    os << tier1_cpi.totalBasecallCount()
       << ':'
       << tier1_cpi.unusedBasecallCount()
       << ':'
       << tier1_cpi.rawPileup().spanningDeletionReadCount
       << ':'
       << tier1_cpi.rawPileup().submappedReadCount;


    std::array<unsigned,N_BASE> tier1_base_counts;
    std::array<unsigned,N_BASE> tier2_base_counts;
    tier1_cpi.cleanedPileup().getBasecallCounts(tier1_base_counts);
    tier2_cpi.cleanedPileup().getBasecallCounts(tier2_base_counts);

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
        const unsigned normalDP(n1_epd.totalBasecallCount());
        const unsigned tumorDP(t1_epd.totalBasecallCount());

        if (dopt.sfilter.is_max_depth())
        {
            if (normalDP > maxChromDepth)
            {
                smod.filters.set(SOMATIC_VARIANT_VCF_FILTERS::HighDepth);
            }
        }

        {
            const unsigned normalFDP(n1_epd.unusedBasecallCount());
            const unsigned tumorFDP(t1_epd.unusedBasecallCount());

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

        // set LowDepth filter if tier 1 tumor or normal depth is below a threshold
        if ((t1_epd.totalBasecallCount() < opt.sfilter.minPassedCallDepth) ||
            (n1_epd.totalBasecallCount() < opt.sfilter.minPassedCallDepth))
        {
            smod.filters.set(SOMATIC_VARIANT_VCF_FILTERS::LowDepth);
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
        os << ";MQ=" << smod.features.get(SOMATIC_SNV_SCORING_FEATURES::RMSMappingQuality);
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

        os << ";ReadPosRankSum=" << smod.features.get(SOMATIC_SNV_SCORING_FEATURES::TumorSampleReadPosRankSum);
        os << ";SNVSB=" << smod.features.get(SOMATIC_SNV_SCORING_FEATURES::TumorSampleStrandBias);

        if (smod.isEVS)
        {
            os << ";" << opt.SomaticEVSVcfInfoTag  << "=" << smod.EVS;
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
    write_vcf_sample_info(dopt,n1_epd,n2_epd,os);

    // tumor sample info:
    os << "\t";
    write_vcf_sample_info(dopt,t1_epd,t2_epd,os);
}

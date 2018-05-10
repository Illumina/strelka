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
///

#include "SomaticIndelVcfWriter.hh"
#include "somaticAlleleUtil.hh"
#include "somatic_call_shared.hh"
#include "somatic_indel_grid.hh"
#include "somatic_indel_scoring_features.hh"
#include "blt_util/blt_exception.hh"
#include "blt_util/io_util.hh"
#include "blt_util/fisher_exact_test.hh"
#include "blt_util/binomial_test.hh"

#include <iomanip>
#include <iostream>



static
void
write_vcf_isri_tiers(
    const AlleleSampleReportInfo& isri1,
    const AlleleSampleReportInfo& isri2,
    const LocalRegionStats& was,
    std::ostream& os)
{
    static const char sep(':');
//  DP:DP2:TAR:TIR:TOR...
    os << isri1.indelLocusDepth
       << sep
       << isri2.indelLocusDepth
       << sep
       << isri1.n_confident_ref_reads+isri1.n_confident_alt_reads << ','
       << isri2.n_confident_ref_reads+isri2.n_confident_alt_reads
       << sep
       << isri1.n_confident_indel_reads << ','
       << isri2.n_confident_indel_reads
       << sep
       << isri1.n_other_reads << ','
       << isri2.n_other_reads;

    const float used(was.regionUsedBasecallCount.avg());
    const float filt(was.regionUnusedBasecallCount.avg());
    const float submap(was.regionSubmappedReadCount.avg());

    const StreamScoper ss(os);
    os << std::fixed << std::setprecision(2);

    os << sep << (used+filt)
       << sep << filt
       << sep << submap
       << sep << calculateBCNoise(was)
       ;
}



static
void
writeSomaticIndelVcfGrid(
    const strelka_options& opt,
    const strelka_deriv_options& dopt,
    const std::string& chromName,
    const pos_t pos,
    const SomaticIndelVcfInfo& siInfo,
    const LocalRegionStats& wasNormal,
    const LocalRegionStats& wasTumor,
    const double maxChromDepth,
    std::ostream& os)
{
    const indel_result_set& rs(siInfo.sindel.rs);

    strelka_shared_modifiers_indel smod;

    // this filter is applied to both EVS and non-EVS models
    /// TODO add a "relative to expected depth" feature to the EVS model so that EVS is a single filter/single ROC model
    if (dopt.sfilter.is_max_depth())
    {
        const unsigned& depth(siInfo.nisri[0].indelLocusDepth);
        if (depth > maxChromDepth)
        {
            smod.filters.set(SOMATIC_VARIANT_VCF_FILTERS::HighDepth);
        }
    }

    // calculate empirical scoring score and features
    calculateScoringFeatures(siInfo, wasNormal, wasTumor, opt, smod);

    const bool is_use_empirical_scoring(opt.isUseSomaticIndelScoring());
    if (!is_use_empirical_scoring)
    {
        smod.EVS = 0;
        smod.isEVS = false;

        // compute all site filters:
        const double normalWinFrac = calculateBCNoise(wasNormal);
        const double tumorWinFrac = calculateBCNoise(wasTumor);

        if ((normalWinFrac >= opt.sfilter.indelMaxWindowFilteredBasecallFrac) ||
            (tumorWinFrac >= opt.sfilter.indelMaxWindowFilteredBasecallFrac))
        {
            smod.filters.set(SOMATIC_VARIANT_VCF_FILTERS::IndelBCNoise);
        }

        if ((rs.ntype != NTYPE::REF) || (rs.from_ntype_qphred < opt.sfilter.sindelQuality_LowerBound))
        {
            smod.filters.set(SOMATIC_VARIANT_VCF_FILTERS::QSI_ref);
        }
    }
    else
    {
        assert(dopt.somaticIndelScoringModel);
        const VariantScoringModelServer& varModel(*dopt.somaticIndelScoringModel);
        updateAlleleEVSScore(varModel, rs, smod);

        if (smod.EVS < varModel.scoreFilterThreshold())
        {
            smod.filters.set(SOMATIC_VARIANT_VCF_FILTERS::LowEVSindel);
        }
    }

    // set LowDepth filter if tier 1 tumor or normal depth is below a threshold
    if ((siInfo.tisri[0].indelLocusDepth < opt.sfilter.minPassedCallDepth) ||
        (siInfo.nisri[0].indelLocusDepth < opt.sfilter.minPassedCallDepth))
    {
        smod.filters.set(SOMATIC_VARIANT_VCF_FILTERS::LowDepth);
    }

    const pos_t output_pos(pos+1);

    static const char sep('\t');
    // CHROM
    os << chromName;

    // POS+
    os << sep << output_pos;

    // ID
    os << sep << ".";

    // REF/ALT
    os << sep << siInfo.vcf_ref_seq
       << sep << siInfo.vcf_indel_seq;

    //QUAL:
    os << sep << ".";

    //FILTER:
    os << sep;
    smod.filters.write(os);

    //INFO
    os << sep
       << "SOMATIC";

    os << ";QSI=" << rs.qphred
       << ";TQSI=" << (siInfo.sindel.sindel_tier+1)
       << ";NT=" << NTYPE::label(rs.ntype)
       << ";QSI_NT=" << rs.from_ntype_qphred
       << ";TQSI_NT=" << (siInfo.sindel.sindel_from_ntype_tier+1)
       << ";SGT="; // << static_cast<SOMATIC_DIGT::index_t>(rs.max_gt);
    DDIGT::write_indel_state(static_cast<DDIGT::index_t>(rs.max_gt),os);

    {
        //MQ and MQ0
        MapqTracker mapqTracker(siInfo.nisri[1].mapqTracker);
        mapqTracker.merge(siInfo.tisri[1].mapqTracker);

        const StreamScoper ss(os);
        os << std::fixed << std::setprecision(2);
        os << ";MQ=" << mapqTracker.getRMS()
           << ";MQ0=" << mapqTracker.zeroCount;

        if (siInfo.indelReportInfo.isRepeatUnit())
        {
            os << ";RU=" << siInfo.indelReportInfo.repeatUnit
               << ";RC=" << siInfo.indelReportInfo.refRepeatCount
               << ";IC=" << siInfo.indelReportInfo.indelRepeatCount;
        }
        os << ";IHP=" << siInfo.indelReportInfo.interruptedHomopolymerLength;

        if (smod.isEVS)
        {
            os << ";" << opt.SomaticEVSVcfInfoTag << "=" << smod.EVS;
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

    // vcf header does not include breakpoint fields, so don't let this be casually turned back on:
    assert (siInfo.indelReportInfo.it != SimplifiedIndelReportType::BREAKPOINT);

    if (rs.is_overlap)
    {
        os << ";OVERLAP";
    }


    //FORMAT
    os << sep << "DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50:BCN50";

    // write normal sample info:
    os << sep;
    write_vcf_isri_tiers(siInfo.nisri[0],siInfo.nisri[1], wasNormal,os);

    // write tumor sample info:
    os << sep;
    write_vcf_isri_tiers(siInfo.tisri[0],siInfo.tisri[1], wasTumor,os);

    os << '\n';
}



void
SomaticIndelVcfWriter::
cacheIndel(
    const pos_t pos,
    const SomaticIndelVcfInfo& siInfo)
{
#if 0
    if (_data.count(pos) != 0)
    {
        std::ostringstream oss;
        oss << "Attempting to cache 2 indels at one site.\n";
        oss << "\texisting indel REF/ALT:\t" << siInfo.iri.ref_seq << " " << siInfo.iri.indel_seq << "\n";
        oss << "\tnew indel REF/ALT:\t" << _data[pos].iri.ref_seq << " " << _data[pos].iri.indel_seq << "\n";
        throw blt_exception(oss.str().c_str());
    }
#endif
    _data[pos].push_back(siInfo);
}



void
SomaticIndelVcfWriter::
addIndelWindowData(
    const std::string& chromName,
    const pos_t pos,
    const LocalRegionStats& wasNormal,
    const LocalRegionStats& wasTumor,
    const double maxChromDepth)
{
    assert(not chromName.empty());
    assert(testPos(pos));
    for (const auto& indelInfo : _data[pos])
    {
        writeSomaticIndelVcfGrid(_opt, _dopt, chromName, pos, indelInfo, wasNormal, wasTumor, maxChromDepth, *_osptr);
    }
    _data.erase(pos);
}

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
///

#include "SomaticIndelVcfWriter.hh"
#include "somatic_call_shared.hh"
#include "somatic_indel_grid.hh"
#include "somatic_indel_scoring_features.hh"
#include "strelka_vcf_locus_info.hh"
#include "blt_util/blt_exception.hh"
#include "blt_util/io_util.hh"
#include "blt_util/qscore.hh"
#include "blt_util/fisher_exact_test.hh"
#include "blt_util/binomial_test.hh"

#include <iomanip>
#include <iostream>
#include <limits>



static
void
write_vcf_isri_tiers(
    const AlleleSampleReportInfo& isri1,
    const AlleleSampleReportInfo& isri2,
    const win_avg_set& was,
    std::ostream& os)
{
    static const char sep(':');
//  DP:DP2:TAR:TIR:TOR...
    os << isri1.tier1Depth
       << sep
       << isri2.tier1Depth
       << sep
       << isri1.n_confident_ref_reads+isri1.n_confident_alt_reads << ','
       << isri2.n_confident_ref_reads+isri2.n_confident_alt_reads
       << sep
       << isri1.n_confident_indel_reads << ','
       << isri2.n_confident_indel_reads
       << sep
       << isri1.n_other_reads << ','
       << isri2.n_other_reads;

    const float used(was.ss_used_win.avg());
    const float filt(was.ss_filt_win.avg());
    const float submap(was.ss_submap_win.avg());

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
    const pos_t pos,
    const SomaticIndelVcfInfo& siInfo,
    const win_avg_set& wasNormal,
    const win_avg_set& wasTumor,
    std::ostream& os)
{
    const indel_result_set& rs(siInfo.sindel.rs);

    // always use high depth filter when enabled
    strelka_shared_modifiers_indel smod;
    if (dopt.sfilter.is_max_depth())
    {
        const unsigned& depth(siInfo.nisri[0].tier1Depth);
        if (depth > dopt.sfilter.max_chrom_depth)
        {
            smod.filters.set(SOMATIC_VARIANT_VCF_FILTERS::HighDepth);
        }
    }

    // calculate empirical scoring score and features
    calculateScoringFeatures(siInfo, wasNormal, wasTumor, opt, dopt, smod);

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
        smod.EVS = varModel.scoreVariant(smod.features.getAll());

        static const int maxEmpiricalVariantScore(60);
        smod.EVS = std::min(int(error_prob_to_phred(smod.EVS)),maxEmpiricalVariantScore);

        smod.isEVS = true;

        if (rs.ntype != NTYPE::REF)
        {
        	smod.EVS = 1.0;
            smod.filters.set(SOMATIC_VARIANT_VCF_FILTERS::Nonref);
        }

        if (smod.EVS < varModel.scoreFilterThreshold())
        {
            smod.filters.set(SOMATIC_VARIANT_VCF_FILTERS::LowEVSindel);
        }
    }

    const pos_t output_pos(pos+1);

    static const char sep('\t');
    // CHROM
    os << opt.bam_seq_name;

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

    if (smod.isEVS)
    {
        const StreamScoper ss(os);
        os << std::fixed << std::setprecision(4);
        os << ";EVS=" << smod.EVS;
    }

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
           << ";MQ0=" << mapqTracker.zeroCount
           ;
    }
    if (siInfo.indelReportInfo.is_repeat_unit())
    {
        os << ";RU=" << siInfo.indelReportInfo.repeat_unit
           << ";RC=" << siInfo.indelReportInfo.ref_repeat_count
           << ";IC=" << siInfo.indelReportInfo.indel_repeat_count;
    }
    os << ";IHP=" << siInfo.indelReportInfo.ihpol;

    if (opt.isReportEVSFeatures)
    {
        const StreamScoper ss(os);
        os << std::setprecision(5);
        os << ";EVSF=";
        smod.features.writeValues(os);
        os << ",";
        smod.dfeatures.writeValues(os);
    }

    if (siInfo.indelReportInfo.it == SimplifiedIndelReportType::BREAKPOINT)
    {
        os << ";SVTYPE=BND";
    }

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
        oss << "ERROR: Attempting to cache 2 indels at one site.\n";
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
    const pos_t pos,
    const win_avg_set& wasNormal,
    const win_avg_set& wasTumor)
{
    assert(testPos(pos));
    for (const auto& indelInfo : _data[pos])
    {
        writeSomaticIndelVcfGrid(_opt, _dopt, pos, indelInfo, wasNormal, wasTumor, *_osptr);
    }
    _data.erase(pos);
}

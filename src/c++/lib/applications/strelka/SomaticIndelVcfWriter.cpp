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

#include "SomaticIndelVcfWriter.hh"
#include "somatic_call_shared.hh"
#include "somatic_indel_grid.hh"
#include "strelka_vcf_locus_info.hh"

#include <iostream>



static
void
write_vcf_isri_tiers(
    const starling_indel_sample_report_info& isri1,
    const starling_indel_sample_report_info& isri2,
    const win_avg_set& was,
    std::ostream& os)
{
    static const char sep(':');
    os << isri1.depth
       << sep
       << isri2.depth
       << sep
       << isri1.n_q30_ref_reads+isri1.n_q30_alt_reads << ','
       << isri2.n_q30_ref_reads+isri2.n_q30_alt_reads
       << sep
       << isri1.n_q30_indel_reads << ','
       << isri2.n_q30_indel_reads
       << sep
       << isri1.n_other_reads << ','
       << isri2.n_other_reads;

    const int used(was.ss_used_win.avg());
    const int filt(was.ss_filt_win.avg());
    const int submap(was.ss_submap_win.avg());
    os << sep << (used+filt)
       << sep << filt
       << sep << submap;
}



static
double
safeFrac(const int num, const int denom)
{
    return ( (denom > 0) ? (num/static_cast<double>(denom)) : 0.);
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
    const somatic_indel_call::result_set& rs(siInfo.sindel.rs);

    strelka_shared_modifiers smod;
    {
        // compute all site filters:
        if (dopt.sfilter.is_max_depth())
        {
            const unsigned& depth(siInfo.nisri[0].depth);
            if (depth > dopt.sfilter.max_depth)
            {
                smod.set_filter(STRELKA_VCF_FILTERS::HighDepth);
            }
        }

        if (siInfo.iri.ref_repeat_count > opt.sfilter.indelMaxRefRepeat)
        {
            smod.set_filter(STRELKA_VCF_FILTERS::Repeat);
        }

        if (siInfo.iri.ihpol > opt.sfilter.indelMaxIntHpolLength)
        {
            smod.set_filter(STRELKA_VCF_FILTERS::iHpol);
        }

        {
            const int normalFilt(wasNormal.ss_filt_win.avg());
            const int normalUsed(wasNormal.ss_used_win.avg());
            const float normalWinFrac(safeFrac(normalFilt,(normalFilt+normalUsed)));

            const int tumorFilt(wasTumor.ss_filt_win.avg());
            const int tumorUsed(wasTumor.ss_used_win.avg());
            const float tumorWinFrac(safeFrac(tumorFilt,(tumorFilt+tumorUsed)));

            if ((normalWinFrac >= opt.sfilter.indelMaxWindowFilteredBasecallFrac) ||
                (tumorWinFrac >= opt.sfilter.indelMaxWindowFilteredBasecallFrac))
            {
                smod.set_filter(STRELKA_VCF_FILTERS::IndelBCNoise);
            }
        }

        if ((rs.ntype != NTYPE::REF) || (rs.sindel_from_ntype_qphred < opt.sfilter.sindelQuality_LowerBound))
        {
            smod.set_filter(STRELKA_VCF_FILTERS::QSI_ref);
        }
    }

    const pos_t output_pos(pos+1);

    // CHROM
    os << opt.bam_seq_name;

    // POS+
    os << '\t' << output_pos;

    // REF/ALT
    os << '\t' << siInfo.iri.vcf_ref_seq
       << '\t' << siInfo.iri.vcf_indel_seq;

    //QUAL:
    os << "\t.";

    //FILTER:
    os << "\t";
    smod.write_filters(os);

    //INFO
    os << '\t'
       << "SOMATIC"
       << ";QSI=" << rs.sindel_qphred
       << ";TQSI=" << (siInfo.sindel.sindel_tier+1)
       << ";NT=" << NTYPE::label(rs.ntype)
       << ";QSI_NT=" << rs.sindel_from_ntype_qphred
       << ";TQSI_NT=" << (siInfo.sindel.sindel_from_ntype_tier+1)
       << ";SGT=" << static_cast<DDIINDEL_GRID::index_t>(rs.max_gt);
    if (siInfo.iri.repeat_unit != "N/A")
    {
        os << ";RU=" << siInfo.iri.repeat_unit
           << ";RC=" << siInfo.iri.ref_repeat_count
           << ";IC=" << siInfo.iri.indel_repeat_count;
    }
    os << ";IHP=" << siInfo.iri.ihpol;
    if ((siInfo.iri.it == INDEL::BP_LEFT) ||
        (siInfo.iri.it == INDEL::BP_RIGHT))
    {
        os << ";SVTYPE=BND";
    }
    if (rs.is_overlap)
    {
        os << ";OVERLAP";
    }


    //FORMAT
    os << '\t' << "DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50";

    // write normal sample info:
    os << '\t';
    write_vcf_isri_tiers(siInfo.nisri[0],siInfo.nisri[1], wasNormal,os);

    // write tumor sample info:
    os << '\t';
    write_vcf_isri_tiers(siInfo.tisri[0],siInfo.tisri[1], wasTumor,os);

    os << '\n';
}


void
SomaticIndelVcfWriter::
cacheIndel(
    const pos_t pos,
    const SomaticIndelVcfInfo& siInfo)
{
    assert(_data.count(pos) == 0);
    _data[pos] = siInfo;
}



void
SomaticIndelVcfWriter::
addIndelWindowData(
    const pos_t pos,
    const win_avg_set& wasNormal,
    const win_avg_set& wasTumor)
{
    assert( _data.count(pos) != 0);
    writeSomaticIndelVcfGrid(_opt, _dopt, pos, _data[pos], wasNormal, wasTumor, *_osptr);

    _data.erase(pos);
}

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
    std::ostream& os)
{
    os << isri1.depth
       << ':'
       << isri2.depth
       << ':'
       << isri1.n_q30_ref_reads+isri1.n_q30_alt_reads << ','
       << isri2.n_q30_ref_reads+isri2.n_q30_alt_reads
       << ':'
       << isri1.n_q30_indel_reads << ','
       << isri2.n_q30_indel_reads
       << ':'
       << isri1.n_other_reads << ','
       << isri2.n_other_reads;
}



static
void
writeSomaticIndelVcfGrid(
    const strelka_options& opt,
    const strelka_deriv_options& dopt,
    const pos_t pos,
    const SomaticIndelVcfInfo& siInfo,
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

        if ((rs.ntype != NTYPE::REF) || (rs.sindel_from_ntype_qphred < opt.sfilter.sindelQuality_LowerBound))
        {
            smod.set_filter(STRELKA_VCF_FILTERS::QSI_ref);
        }
    }

    // CHROM
    os << opt.bam_seq_name;

    // POS
    os << '\t' << pos;

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
    os << '\t' << "DP:DP2:TAR:TIR:TOR";

    // write normal sample info:
    os << '\t';
    write_vcf_isri_tiers(siInfo.nisri[0],siInfo.nisri[1],os);

    // write tumor sample info:
    os << '\t';
    write_vcf_isri_tiers(siInfo.tisri[0],siInfo.tisri[1],os);

    os << '\n';
}


void
SomaticIndelVcfWriter::
queueIndel(
    const pos_t pos,
    const SomaticIndelVcfInfo& siInfo)
{
    writeSomaticIndelVcfGrid(_opt, _dopt, pos, siInfo, *_osptr);
}

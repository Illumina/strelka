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

#include "denovo_indel_call_vcf.hh"
#include "pedicure_vcf_locus_info.hh"
#include "blt_util/blt_exception.hh"
#include "blt_util/io_util.hh"

#include <iomanip>
#include <iostream>



static
void
write_vcf_sample_info(
    const starling_indel_sample_report_info& isri1,
    const starling_indel_sample_report_info& /*isri2*/,
    const denovo_indel_call& dinc,
    const int& id,
    std::ostream& os)
{

//	GT:GQ:GQX:DP:DP2:AD:PL
    int total = isri1.n_q30_ref_reads + isri1.n_q30_indel_reads;

    static const char sep(':');
    if (dinc.gtstring.size()>0 && total>0)
    {
        os << dinc.gtstring.at(id)
           << sep
           << dinc.gq.at(id)
           << sep
           << dinc.gqx.at(id)
           << sep;
    }
    else
    {
        os << "./." << sep << "." << sep << "." << sep;
    }

    os << isri1.tier1Depth
       << sep
//	   n_q30_alt_reads + n_q30_indel_reads + n_q30_ref_reads;
       << isri1.n_q30_ref_reads << ','
       << isri1.n_q30_indel_reads;// + isri1.n_q30_alt_reads;
//       << isri1.n_q30_ref_reads+isri2.n_q30_ref_reads << ','
//        << isri1.n_q30_alt_reads+isri2.n_q30_alt_reads
//       << sep
//       << isri1.n_q30_indel_reads << ','
//       << isri2.n_q30_indel_reads
//       << sep
//       << isri1.n_other_reads << ','
//       << isri2.n_other_reads;
}


#if 0
static
double
safeFrac(const int num, const int denom)
{
    return ( (denom > 0) ? (num/static_cast<double>(denom)) : 0.);
}
#endif



void
denovo_indel_call_vcf(
    const pedicure_options& opt,
    const pedicure_deriv_options& dopt,
    const SampleInfoManager& sinfo,
    const denovo_indel_call& dinc,
    const starling_indel_report_info& iri,
    const std::vector<isriTiers_t>& isri,
    std::ostream& os)
{
    const denovo_indel_call::result_set& rs(dinc.rs);

    pedicure_shared_modifiers smod;
    {
        // compute all site filters:
        if (dopt.dfilter.is_max_depth())
        {
            using namespace PEDICURE_SAMPLETYPE;
            const unsigned probandIndex(sinfo.getTypeIndexList(PROBAND)[0]);

            const unsigned& depth(isri[probandIndex][PEDICURE_TIERS::TIER1].tier1Depth);
            if (depth > dopt.dfilter.max_depth)
            {
                smod.set_filter(PEDICURE_VCF_FILTERS::HighDepth);
            }
        }

//        if (rs.is_overlap)
//        {
//            smod.set_filter(PEDICURE_VCF_FILTERS::OverlapConflict);
//        }

        for (unsigned sampleIndex(0); sampleIndex<sinfo.size(); sampleIndex++)
        {
            if (dinc.gqx[sampleIndex] < opt.dfilter.sindelQuality_LowerBound)
            {
                smod.set_filter(PEDICURE_VCF_FILTERS::LowGQX);
                break;
            }
        }


        if (rs.dindel_qphred < opt.dfilter.dindel_qual_lowerbound)
        {
//            smod.set_filter(PEDICURE_VCF_FILTERS::QDI);
        }

        if (iri.ref_repeat_count > opt.dfilter.indelMaxRefRepeat)
        {
            smod.set_filter(PEDICURE_VCF_FILTERS::Repeat);
        }

        if (iri.ihpol > opt.dfilter.indelMaxIntHpolLength)
        {
            smod.set_filter(PEDICURE_VCF_FILTERS::iHpol);
        }


//        {
//            const int normalFilt(wasNormal.ss_filt_win.avg());
//            const int normalUsed(wasNormal.ss_used_win.avg());
//            const float normalWinFrac(safeFrac(normalFilt,(normalFilt+normalUsed)));
//
//            const int tumorFilt(wasTumor.ss_filt_win.avg());
//            const int tumorUsed(wasTumor.ss_used_win.avg());
//            const float tumorWinFrac(safeFrac(tumorFilt,(tumorFilt+tumorUsed)));
//
//            if ((normalWinFrac >= opt.sfilter.indelMaxWindowFilteredBasecallFrac) ||
//                (tumorWinFrac >= opt.sfilter.indelMaxWindowFilteredBasecallFrac))
//            {
//                smod.set_filter(STRELKA_VCF_FILTERS::IndelBCNoise);
//            }
//        }

//        if ((rs.ntype != NTYPE::REF) || (rs.sindel_from_ntype_qphred < opt.sfilter.sindelQuality_LowerBound))
//        {
//            smod.set_filter(STRELKA_VCF_FILTERS::QSI_ref);
//        }

    }


    static const char sep('\t');

    // REF/ALT
    os << sep << iri.vcf_ref_seq
       << sep << iri.vcf_indel_seq;

    //QUAL:
    os << sep << ".";

    //FILTER:
    os << sep;
    smod.write_filters(os);

    //INFO
    unsigned totalDepth(0);
    for (const auto& sampleIsri : isri)
    {
        totalDepth += sampleIsri[PEDICURE_TIERS::TIER1].tier1Depth;
    }

    os << sep
       << "DP="<< totalDepth
       << ";DQ=" << rs.dindel_qphred
//       << ";TQDI=" << (dinc.dindel_tier+1)
       ;

    if (iri.is_repeat_unit())
    {
        os << ";RU=" << iri.repeat_unit
           << ";RC=" << iri.ref_repeat_count
           << ";IC=" << iri.indel_repeat_count;
    }
    os << ";IHP=" << iri.ihpol;
//    if ((iri.it == INDEL::BP_LEFT) ||
//        (iri.it == INDEL::BP_RIGHT))
//    {
//        os << ";SVTYPE=BND";
//    }
    if (rs.is_overlap)
    {
//        os << ";OVERLAP";
    }

    //FORMAT
    os << sep << "GT:GQ:GQX:DPI:AD";

    // write sample info:
    int id = 0;
    for (const auto& sampleIsri : isri)
    {
        os << sep;
        write_vcf_sample_info(sampleIsri[PEDICURE_TIERS::TIER1],sampleIsri[PEDICURE_TIERS::TIER2],dinc, id, os);
        id++;
    }
}

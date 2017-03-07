// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2017 Illumina, Inc.
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
    const AlleleSampleReportInfo& indelSampleReportInfo1,
    const AlleleSampleReportInfo& /*indelSampleReportInfo2*/,
    const denovo_indel_call& dinc,
    const int& id,
    std::ostream& os)
{

//	GT:GQ:GQX:DP:DP2:AD:PL
    int total = indelSampleReportInfo1.n_confident_ref_reads + indelSampleReportInfo1.n_confident_indel_reads;

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

    os << indelSampleReportInfo1.tier1Depth
       << sep
       << indelSampleReportInfo1.n_confident_ref_reads << ','
       << indelSampleReportInfo1.n_confident_indel_reads;
}



void
denovo_indel_call_vcf(
    const pedicure_options& opt,
    const pedicure_deriv_options& dopt,
    const SampleInfoManager& sinfo,
    const double maxChromDepth,
    const denovo_indel_call& dinc,
    const AlleleReportInfo& indelReportInfo,
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
            if (depth > maxChromDepth)
            {
                smod.set_filter(PEDICURE_VCF_FILTERS::HighDepth);
            }
        }


        for (unsigned sampleIndex(0); sampleIndex<sinfo.size(); sampleIndex++)
        {
            if (dinc.gqx[sampleIndex] < opt.dfilter.sindelQuality_LowerBound)
            {
                smod.set_filter(PEDICURE_VCF_FILTERS::LowGQX);
                break;
            }
        }

#if 0
        if (rs.dindel_qphred < opt.dfilter.dindel_qual_lowerbound)
        {
            smod.set_filter(PEDICURE_VCF_FILTERS::QDI);
        }
#endif

        if (indelReportInfo.ref_repeat_count > opt.dfilter.indelMaxRefRepeat)
        {
            smod.set_filter(PEDICURE_VCF_FILTERS::Repeat);
        }

        if (indelReportInfo.ihpol > opt.dfilter.indelMaxIntHpolLength)
        {
            smod.set_filter(PEDICURE_VCF_FILTERS::iHpol);
        }
    }


    static const char sep('\t');

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

    if (indelReportInfo.is_repeat_unit())
    {
        os << ";RU=" << indelReportInfo.repeat_unit
           << ";RC=" << indelReportInfo.ref_repeat_count
           << ";IC=" << indelReportInfo.indel_repeat_count;
    }
    os << ";IHP=" << indelReportInfo.ihpol;

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

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

#include "denovo_snv_call_vcf.hh"
#include "pedicure_vcf_locus_info.hh"
#include "blt_util/io_util.hh"

#include <array>
#include <iomanip>
#include <iostream>

#include "blt_util/log.hh"

static
void
write_vcf_sample_info(
    const blt_options& opt,
    const CleanedPileup& tier1_cpi,
    const CleanedPileup& tier2_cpi,
    const denovo_snv_call& dsc,
    int sampleIndex,
    std::ostream& os)
{
    //DP:FDP:SDP:SUBDP:AU:CU:GU:TU
//	log_os << dsc.Sampleplhoods[sampleIndex][0]
//			<< " " << dsc.Sampleplhoods[sampleIndex][1]
//			<< " " << dsc.Sampleplhoods[sampleIndex][2] <<	std::endl;
	if (dsc.gts[sampleIndex]==0)
    	os << "0/0";
	if (dsc.gts[sampleIndex]==1)
		os << "0/1";
	if (dsc.gts[sampleIndex]==2)
		os << "1/1";
	os <<':'
	   << "1" //GQ
	   <<':'
	   << dsc.gqx[sampleIndex]  //GQX
       <<':'
	   << tier1_cpi.n_calls()
       << ':'
       << tier1_cpi.n_unused_calls()
       << ':'
       << "0,23"
       << ':'
       << "1,2";
//       << tier1_cpi.rawPileup().n_spandel
//       << ':'
//       << tier1_cpi.rawPileup().n_submapped;

    std::array<unsigned,N_BASE> tier1_base_counts;
    std::array<unsigned,N_BASE> tier2_base_counts;
    tier1_cpi.cleanedPileup().get_known_counts(tier1_base_counts,opt.used_allele_count_min_qscore);
    tier2_cpi.cleanedPileup().get_known_counts(tier2_base_counts,opt.used_allele_count_min_qscore);
}

void
denovo_snv_call_vcf(
    const pedicure_options& opt,
    const pedicure_deriv_options& dopt,
    const SampleInfoManager& sinfo,
    const cpiPtrTiers_t& pileups,
    const denovo_snv_call& dsc,
    std::ostream& os)
{
    using namespace PEDICURE_SAMPLETYPE;

    const unsigned probandIndex(sinfo.getTypeIndexList(PROBAND)[0]);
    const CleanedPileup& probandCpi(*pileups[PEDICURE_TIERS::TIER1][probandIndex]);

    const denovo_snv_call::result_set& rs(dsc.rs);

    pedicure_shared_modifiers smod;

    {
        // compute all site filters:
        const unsigned probandDP(probandCpi.n_calls());

        if (dopt.dfilter.is_max_depth())
        {
            if (probandDP > dopt.dfilter.max_depth)
            {
                smod.set_filter(PEDICURE_VCF_FILTERS::HighDepth);
            }
        }

        if (rs.dsnv_qphred < opt.dfilter.dsnv_qual_lowerbound)
        {
//            smod.set_filter(PEDICURE_VCF_FILTERS::QDS);
        }

    }


    //REF:
    os << '\t' << probandCpi.rawPileup().get_ref_base()
       //ALT:
       << "\t"
	   << dsc.alt_str;
    //QUAL:
    os << "\t.";

    //FILTER:
    os << "\t";
    smod.write_filters(os);

    //INFO:
    os << "\t";

    {
        const StreamScoper ss(os);
        os << std::fixed << std::setprecision(2);

        // m_mapq includes all calls, even from reads below the mapq threshold:
        unsigned n_mapq(0);
        unsigned n_mapq0(0);
        for (unsigned sampleIndex(0); sampleIndex<sinfo.size(); ++sampleIndex)
        {
            const snp_pos_info& pi(pileups[PEDICURE_TIERS::TIER1][sampleIndex]->rawPileup());
            n_mapq += pi.n_mapq;
            n_mapq0 += pi.n_mapq0;
        }
        os << "DP=" << n_mapq;
        os << ";MQ0=" << n_mapq0;
        os << ";QDS=" << rs.dsnv_qphred;
//           << ";TQSI=" << (dsc.dsnv_tier+1);

    }

    //FORMAT:
    os << '\t'
       << "GT:GQ:GQX:DP:FDP:AD:PL";

    for (unsigned sampleIndex(0); sampleIndex<sinfo.size(); sampleIndex++)
    {
        const CleanedPileup& cpi1(*pileups[PEDICURE_TIERS::TIER1][sampleIndex]);
        const CleanedPileup& cpi2(*pileups[PEDICURE_TIERS::TIER2][sampleIndex]);
        os << "\t";
        write_vcf_sample_info(opt,cpi1,cpi2,dsc,sampleIndex,os);
    }
}

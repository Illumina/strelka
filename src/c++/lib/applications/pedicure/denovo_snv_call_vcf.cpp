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

#include "denovo_snv_call_vcf.hh"
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
    const CleanedPileup& /*tier2_cpi*/,
    const denovo_snv_call& dsc,
    int sampleIndex,
    std::ostream& os)
{

	std::array<unsigned,N_BASE> tier1_base_counts;
	tier1_cpi.cleanedPileup().get_known_counts(tier1_base_counts,opt.used_allele_count_min_qscore);


	//os << dsc.gts_chrom[sampleIndex][0]
	//   << "/"
	//   << dsc.gts_chrom[sampleIndex][1];
    os << dsc.gtstring[sampleIndex]; 
    
	os <<':'
	   << "1" //GQ
	   <<':'
	   << dsc.gqx[sampleIndex]  //GQX
       <<':'
	   << (tier1_cpi.n_calls()-tier1_cpi.n_unused_calls())
       << ':'
       << tier1_cpi.n_unused_calls()
       << ':'
       << tier1_base_counts[dsc.ref_gt];
	   for (unsigned i=0; i<dsc.alts.size(); i++)
		   os << "," << tier1_base_counts[dsc.alts[i]];

	   // PL field
	   //os << ':'
       //<< dsc.get_pl(sampleIndex);
       os << ':' << dsc.Sampleplhoods[sampleIndex][0] << "," << dsc.Sampleplhoods[sampleIndex][1] << "," << dsc.Sampleplhoods[sampleIndex][2];
	   for(unsigned i=3;i<dsc.Sampleplhoods[sampleIndex].size();++i){
		os << "," << dsc.Sampleplhoods[sampleIndex][i];
	   }
	
}

void
denovo_snv_call_vcf(
    const pedicure_options& opt,
    const pedicure_deriv_options& dopt,
    const SampleInfoManager& sinfo,
    const cpiPtrTiers_t& pileups,
    denovo_snv_call& dsc,
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
    //dsc.get_alt();


    //REF:
    os << '\t' << probandCpi.rawPileup().get_ref_base();
       //ALT:
//       << "\t"
//	   << dsc.alt_str;
	os << "\t" << id_to_base(dsc.alts[0]);
    for(unsigned i=1;i<dsc.alts.size(); ++i){ os << "," << id_to_base(dsc.alts[i]); }
    
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

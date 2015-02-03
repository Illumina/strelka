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
#include "inovo_vcf_locus_info.hh"

#include <array>
#include <iomanip>
#include <iostream>



static
void
write_vcf_sample_info(
    const blt_options& opt,
    const CleanedPileup& cpi,
    std::ostream& os)
{
    //DP:FDP:SDP:SUBDP:AU:CU:GU:TU
    os << cpi.n_calls()
       << ':'
       << cpi.n_unused_calls()
       << ':'
       << cpi.rawPileup().n_spandel
       << ':'
       << cpi.rawPileup().n_submapped;

    std::array<unsigned,N_BASE> base_counts;
    cpi.cleanedPileup().get_known_counts(base_counts,opt.used_allele_count_min_qscore);

    for (unsigned b(0); b<N_BASE; ++b)
    {
        os << ':'
           << base_counts[b];
    }
}



void
denovo_snv_caller_vcf(
    const inovo_options& opt,
    const inovo_deriv_options& dopt,
    const SampleInfoManager& sinfo,
    const std::vector<const CleanedPileup*>& pileups,
    const denovo_snv_call& dsc,
    std::ostream& os)
{
    using namespace INOVO_SAMPLETYPE;

    const unsigned probandIndex(sinfo.getTypeIndexList(PROBAND)[0]);
    const CleanedPileup& probandCpi(*pileups[probandIndex]);

    const result_set& rs(dsc.rs);

    inovo_shared_modifiers smod;

    {
        // compute all site filters:
        const unsigned probandDP(probandCpi.n_calls());

        if (dopt.dfilter.is_max_depth())
        {
            if (probandDP > dopt.dfilter.max_depth)
            {
                smod.set_filter(INOVO_VCF_FILTERS::HighDepth);
            }
        }
    }

    //REF:
    os << '\t' << probandCpi.rawPileup().get_ref_base()
       //ALT:
       << "\t" << '.';
//    DDIGT_SGRID::write_alt_alleles(static_cast<DDIGT_SGRID::index_t>(rs.max_gt),
 //                                  dsc.ref_gt,os);
    //QUAL:
    os << "\t.";

    //FILTER:
    os << "\t";
    smod.write_filters(os);

    //INFO:
    os << '\t.';

    {
        const StreamScoper ss(os);
        os << std::fixed << std::setprecision(2);

        // m_mapq includes all calls, even from reads below the mapq threshold:
        unsigned n_mapq(0);
        unsigned n_mapq0(0);
        for (unsigned sampleIndex(0); sampleIndex<sinfo.size(); ++sampleIndex)
        {
            const snp_pos_info& pi(pileups[sampleIndex]->rawPileup());
            n_mapq += pi.n_mapq;
            n_mapq0 += pi.n_mapq0;
        }
        os << ";DP=" << n_mapq;
        os << ";MQ0=" << n_mapq0;
    }

    //FORMAT:
    os << '\t'
       << "DP:FDP:SDP:SUBDP:AU:CU:GU:TU";

    for (unsigned sampleIndex(0);sampleIndex<sinfo.size();sampleIndex++)
    {
        const CleanedPileup& cpi(*pileups[sampleIndex]);
        os << "\t";
        write_vcf_sample_info(opt,cpi,os);
    }
}

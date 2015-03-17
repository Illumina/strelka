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

#include "denovo_snv_caller.hh"
#include "denovo_snv_grid_states.hh"
#include "strelka_common/position_snp_call_grid_lhood_cached.hh"

#include <array>

#if 0

static
bool
checkME(
    const uint8_t a1,
    const uint8_t a2,
    const std::array<unsigned,N_BASE>& p1,
    const std::array<unsigned,N_BASE>& p2)
{
    static const unsigned minCount(2);
    return ((p1[a1]>=minCount) && (p2[a2]>=minCount));
}

#endif


typedef std::array<blt_float_t,DIGT_DGRID::SIZE> dsnv_state_t;



void
get_denovo_snv_call(
    const inovo_options& opt,
    const SampleInfoManager& sinfo,
    const cpiPtrTiers_t& pileups,
    denovo_snv_call& dsc)
{
    using namespace INOVO_SAMPLETYPE;

    const unsigned probandIndex(sinfo.getTypeIndexList(PROBAND)[0]);
    const std::vector<unsigned>& parentIndex(sinfo.getTypeIndexList(PARENT));

    std::vector<unsigned> allIndex(parentIndex);
    allIndex.push_back(probandIndex);

    // escape in case of low sample depth:
    // depth must be at least minDepth in all of proband and parents
    {
        static const unsigned minDepth(10);
        for (const auto sampleIndex : allIndex)
        {
            const CleanedPileup& cpi(*pileups[INOVO_TIERS::TIER1][sampleIndex]);
            if (cpi.cleanedPileup().calls.size() < minDepth)
            {
                return;
            }
        }
    }

    // setup ref GT
    const CleanedPileup& probandCpi(*pileups[INOVO_TIERS::TIER1][probandIndex]);
    const char refBase(probandCpi.cleanedPileup().get_ref_base());
    if (refBase=='N')
    {
        return;
    }
    dsc.ref_gt=base_to_id(refBase);

    const unsigned sampleSize(sinfo.size());
    std::vector<dsnv_state_t> sampleLhood(sampleSize);


    std::array<denovo_snv_call::result_set,INOVO_TIERS::SIZE> tier_rs;
    for (unsigned tierIndex(0); tierIndex<INOVO_TIERS::SIZE; ++tierIndex)
    {
        const bool is_include_tier2(tierIndex==1);

        if (is_include_tier2)
        {
            if (! opt.tier2.is_tier2()) continue;
            if (tier_rs[0].snv_qphred==0)
            {
                if (! dsc.is_forced_output)   // if forced output then there's still a point to computing tier2
                {
                    tier_rs[1].snv_qphred=0;
                    continue;
                }
            }
        }

        for (unsigned sampleIndex(0);sampleIndex<sampleSize;++sampleIndex)
        {
            // get likelihoods for each sample
            const CleanedPileup& cpi(*pileups[tierIndex][sampleIndex]);
            const snp_pos_info& pi(cpi.cleanedPileup());
            blt_float_t* lhood(sampleLhood[sampleIndex].data());
            get_diploid_gt_lhood_cached(opt, pi, lhood);
            get_diploid_het_grid_lhood_cached(pi, DIGT_DGRID::HET_RES, lhood);

        }

        //calculate_result_set();
    }



    if (! dsc.is_forced_output)
    {
        for (const auto& val : tier_rs)
        {
            if (val.snv_qphred==0) return;
        }
    }

    dsc.snv_tier=0;
    if (opt.tier2.is_tier2())
    {
        if (tier_rs[0].snv_qphred > tier_rs[1].snv_qphred)
        {
            dsc.snv_tier=1;
        }
    }

    dsc.rs=tier_rs[dsc.snv_tier];


#if 0

    // top two alleles;
    uint8_t max1(0), max2(0);
    {
        const CleanedPileup& probandCpi(*pileups[probandIndex]);

        std::array<unsigned,N_BASE> base_count;
        probandCpi.cleanedPileup().get_known_counts(base_count,opt.min_qscore);

        unsigned total(0);
        for (const unsigned count : base_count)
        {
            total += count;
        }
        const unsigned min_allele_count(total*.25);

        for (unsigned i(0);i<N_BASE;++i)
        {
            if(base_count[i] > base_count[max1]) max1=i;
        }

        bool is_max2(false);
        for (unsigned i(0);i<N_BASE;++i)
        {
            if (i==max1) continue;
            if (is_max2 && (base_count[i] <= base_count[max2])) continue;
            if (base_count[i]<min_allele_count) continue;
            max2=i;
            is_max2=true;
        }

        if (! is_max2) max2 = max1;
    }

    // the nature of the hack count model requires two parents:
    if (parentIndex.size() < 2) return;

    //
    {
        const CleanedPileup& p1Cpi(*pileups[parentIndex[0]]);
        const CleanedPileup& p2Cpi(*pileups[parentIndex[1]]);

        std::array<unsigned,N_BASE> p1Count;
        p1Cpi.cleanedPileup().get_known_counts(p1Count,opt.min_qscore);

        std::array<unsigned,N_BASE> p2Count;
        p2Cpi.cleanedPileup().get_known_counts(p2Count,opt.min_qscore);

        // scenario 1: max1 -> parent1, max2 -> parent2
        bool isValid(checkME(max1,max2,p1Count,p2Count));
        if (! isValid)
        {
            // scenario 2: max1 -> parent2, max2 -> parent1
            isValid = checkME(max1,max2,p2Count,p1Count);
        }
        if (isValid)
        {
            return;
        }
    }

    dsc.rs.isCall = true;
#endif
}


#if 0
// strawman model treats normal and tumor as independent, so
// calculate separate lhoods:
blt_float_t normal_lhood[DIGT_SGRID::SIZE];
blt_float_t tumor_lhood[DIGT_SGRID::SIZE];

const bool is_tier2(NULL != normal_epi_t2_ptr);

static const unsigned n_tier(2);
result_set tier_rs[n_tier];
for (unsigned i(0); i<n_tier; ++i)
{
    const bool is_include_tier2(i==1);
    if (is_include_tier2)
    {
        if (! is_tier2) continue;
        if (tier_rs[0].snv_qphred==0)
        {
            tier_rs[1] = tier_rs[0];
            continue;
        }
    }

    // get likelihood of each genotype
    //
    static constexpr bool is_normal_het_bias(false);
    static constexpr blt_float_t normal_het_bias(0.0);
    static constexpr bool is_tumor_het_bias(false);
    static constexpr blt_float_t tumor_het_bias(0.0);

    const extended_pos_info& nepi(is_include_tier2 ? *normal_epi_t2_ptr : normal_epi );
    const extended_pos_info& tepi(is_include_tier2 ? *tumor_epi_t2_ptr : tumor_epi );
    get_diploid_gt_lhood_cached(_opt,nepi.pi,is_normal_het_bias,normal_het_bias,normal_lhood);
    get_diploid_gt_lhood_cached(_opt,tepi.pi,is_tumor_het_bias,tumor_het_bias,tumor_lhood);

    get_diploid_het_grid_lhood_cached(nepi.pi, DIGT_SGRID::HET_RES, normal_lhood+DIGT::SIZE);
    get_diploid_het_grid_lhood_cached(tepi.pi, DIGT_SGRID::HET_RES, tumor_lhood+DIGT::SIZE);

    get_diploid_strand_grid_lhood_spi(nepi.pi,sgt.ref_gt,normal_lhood+DIGT_SGRID::PRESTRAND_SIZE);
    get_diploid_strand_grid_lhood_spi(tepi.pi,sgt.ref_gt,tumor_lhood+DIGT_SGRID::PRESTRAND_SIZE);

    // genomic site results:
    calculate_result_set_grid(isComputeNonSomatic,
                              normal_lhood,
                              tumor_lhood,
                              get_prior_set(sgt.ref_gt),
                              _ln_som_match,_ln_som_mismatch,
                              sgt.ref_gt,
                              sgt.is_forced_output,
                              tier_rs[i]);


}

if (! (sgt.is_forced_output || isComputeNonSomatic))
{
    if ((tier_rs[0].snv_qphred==0) ||
        (is_tier2 && (tier_rs[1].snv_qphred==0))) return;
}

sgt.snv_tier=0;
sgt.snv_from_ntype_tier=0;
if (is_tier2)
{
    if (tier_rs[0].snv_qphred > tier_rs[1].snv_qphred)
    {
        sgt.snv_tier=1;
    }

    if (tier_rs[0].snv_from_ntype_qphred > tier_rs[1].snv_from_ntype_qphred)
    {
        sgt.snv_from_ntype_tier=1;
    }
}

sgt.rs=tier_rs[sgt.snv_from_ntype_tier];

if (is_tier2 && (tier_rs[0].ntype != tier_rs[1].ntype))
{
    // catch NTYPE conflict states:
    sgt.rs.ntype = NTYPE::CONFLICT;
    sgt.rs.snv_from_ntype_qphred = 0;
}
else
{
    // classify NTYPE:
    //

    // convert diploid genotype into more limited ntype set:
    //
    if       (sgt.rs.ntype==sgt.ref_gt)
    {
        sgt.rs.ntype=NTYPE::REF;
    }
    else if (DIGT::is_het(sgt.rs.ntype))
    {
        sgt.rs.ntype=NTYPE::HET;
    }
    else
    {
        sgt.rs.ntype=NTYPE::HOM;
    }
}

sgt.rs.snv_qphred = tier_rs[sgt.snv_tier].snv_qphred;

/// somatic gVCF, always use tier1 to keep things simple:
sgt.rs.nonsomatic_qphred = tier_rs[0].nonsomatic_qphred;
}

#endif

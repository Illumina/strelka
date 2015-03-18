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
#include "blt_util/math_util.hh"
#include "blt_util/prob_util.hh"
#include "blt_util/qscore.hh"
#include "strelka_common/position_snp_call_grid_lhood_cached.hh"

#include <array>
#include <iterator>


//#define DENOVO_SNV_DEBUG

#ifdef DENOVO_SNV_DEBUG
#include "blt_util/log.hh"
#endif



typedef std::array<blt_float_t,DIGT_DGRID::SIZE> dsnv_state_t;



namespace TRANSMISSION_STATE
{
    // "ERROR" represents a de-novo event that is incredibly unlikely (multiple events)
    //  -- we could also put it in the denovo state and just use the denovo prior
    // squared to get the same result -- then the dominant term would actually be the
    // probably of an erroneous copy number observation in the sample instead.
    enum index_t
    {
        INHERITED,
        DENOVO,
        ERROR,
        SIZE
    };

#ifdef DENOVO_SNV_DEBUG
    static
    const char*
    getLabel(
        const index_t idx)
    {
        switch (idx)
        {
        case INHERITED: return "INHERITED";
        case DENOVO: return "DENOVO";
        case ERROR: return "ERROR";
        default:
            assert(false && "Unknown transmission state");
            return nullptr;
        }
    }
#endif

    // temporary fixed priors:
    static
    double
    getPrior(
        const index_t idx)
    {
        // as currently defined background exp is I: 15/27 E: 2/27 D: 10/27 -- compared to drate this doesn't matter
        static const double lndrate(std::log(1e-6));
        static const double noiserate(std::log(1e-7));
        switch (idx)
        {
        case INHERITED: return 0.;
        case DENOVO: return lndrate;
        case ERROR: return noiserate;
        default:
            assert(false && "Undefined inheritance state");
        }
    }

    static const unsigned alleleCount(2);
    typedef std::array<uint8_t,alleleCount> alleleSet_t;

    static
    unsigned
    getTransmissionErrorCount(
        const alleleSet_t& pA,
        const alleleSet_t& pB,
        const alleleSet_t& c)
    {
        unsigned val(0);
        if ((c[0] != pA[0]) && (c[0] != pA[1])) val += 1;
        if ((c[1] != pB[0]) && (c[1] != pB[1])) val += 1;
        return val;
    }

    static
    index_t
    get_state(
        const unsigned parent0GT,
        const unsigned parent1GT,
        const unsigned probandGT)
    {
        alleleSet_t parent0Alleles;
        alleleSet_t parent1Alleles;
        alleleSet_t probandAlleles;
        for (unsigned alleleIndex(0); alleleIndex<alleleCount; ++alleleIndex)
        {
            parent0Alleles[alleleIndex] = DIGT::get_allele(parent0GT,alleleIndex);
            parent1Alleles[alleleIndex] = DIGT::get_allele(parent1GT,alleleIndex);
            probandAlleles[alleleIndex] = DIGT::get_allele(probandGT,alleleIndex);
        }
        const unsigned errorCount(std::min(getTransmissionErrorCount(parent0Alleles,parent1Alleles,probandAlleles),getTransmissionErrorCount(parent1Alleles,parent0Alleles,probandAlleles)));
        switch (errorCount)
        {
        case 0: return INHERITED;
        case 1: return DENOVO;
        case 2: return ERROR;
        default:
            assert(false && "Unexpected error count value");
            return ERROR;
        }
    }
}



template <typename Iter>
typename std::iterator_traits<Iter>::difference_type
max_element_index(Iter b, Iter e)
{
    return std::distance(b, std::max_element(b,e));
}



/// Given the likelihood, go through the final computations to get the
/// posterior and derived values.
///
static
void
calculate_result_set_grid(
    const SampleInfoManager& sinfo,
    const std::vector<dsnv_state_t>& sampleLhood,
    denovo_snv_call::result_set& rs)
{
    using namespace INOVO_SAMPLETYPE;

    const unsigned probandIndex(sinfo.getTypeIndexList(PROBAND)[0]);
    const std::vector<unsigned>& parentIndex(sinfo.getTypeIndexList(PARENT));

    static const double lnzero(-std::numeric_limits<double>::infinity());
    std::array<double,TRANSMISSION_STATE::SIZE> stateLhood;
    std::fill(stateLhood.begin(),stateLhood.end(),lnzero);

    // partial brute force enumeration of all parent-child genotypes:
    //  find max_gt from parents and children, translate this into a candidate allele pool
    //   enumerate all genotypes in all samples from the candidate allele pool only.
    std::vector<unsigned> trioIndex(parentIndex);
    trioIndex.push_back(probandIndex);

    std::array<bool,N_BASE> max_alleles;
    std::fill(max_alleles.begin(),max_alleles.end(),false);
    for (const auto sampleIndex : trioIndex)
    {
        const auto& lhood(sampleLhood[sampleIndex]);
        const unsigned maxGT(max_element_index(lhood.begin(),lhood.begin()+DIGT::SIZE));
        for (unsigned chromCopyIndex(0); chromCopyIndex<2; ++chromCopyIndex)
        {
            max_alleles[DIGT::get_allele(maxGT,chromCopyIndex)] = true;
        }
    }

    auto isFilterGT = [&](const unsigned gt)
        {
            return (! (max_alleles[DIGT::get_allele(gt,0)] && max_alleles[DIGT::get_allele(gt,1)]));
        };

    for (unsigned parent0GT(0); parent0GT<DIGT::SIZE; ++parent0GT)
    {
        if (isFilterGT(parent0GT)) continue;
        for (unsigned parent1GT(0); parent1GT<DIGT::SIZE; ++parent1GT)
        {
            if (isFilterGT(parent1GT)) continue;
            for (unsigned probandGT(0); probandGT<DIGT::SIZE; ++probandGT)
            {
                if (isFilterGT(probandGT)) continue;
                const double pedigreeLhood = sampleLhood[parentIndex[0]][parent0GT] + sampleLhood[parentIndex[1]][parent1GT] + sampleLhood[probandIndex][probandGT];
                const TRANSMISSION_STATE::index_t tran(TRANSMISSION_STATE::get_state(parent0GT,parent1GT,probandGT));
#ifdef DENOVO_SNV_DEBUG
                {
                    using namespace TRANSMISSION_STATE;
                    log_os << "p0gt/p1/c: "
                            << DIGT::get_gt_label(parent0GT) << " "
                            << DIGT::get_gt_label(parent1GT) << " "
                            << DIGT::get_gt_label(probandGT) << " trans_state: " << getLabel(tran) << " lhood: " << pedigreeLhood << "\n";
                }
#endif
                stateLhood[tran] = log_sum(stateLhood[tran],pedigreeLhood);
            }
        }
    }

    // apart from the regular genotype analysis, we go through a bunch of noise states and
    // dump these into the "error' transmission state
    //
    // these are all non-standard allele frequencies shared among all samples --
    // 99% of the time this is meant to catch low-frequency alt noise shared in all three
    // samples (if all were sequenced to a very high depth) but spuriously more prevalent
    // in the proband due to low-depth sampling issues:
    static const unsigned ratioCount(DIGT_DGRID::HET_RES*2);
    for (unsigned hetIndex(0); hetIndex<(DIGT::HET_SIZE); hetIndex++)
    {
        const unsigned hetGT(hetIndex+N_BASE);
        if (isFilterGT(hetGT)) continue;
        for (unsigned ratioIndex(0); ratioIndex<ratioCount; ++ratioIndex)
        {
            const unsigned noiseState(DIGT::SIZE+(ratioIndex*DIGT::HET_SIZE)+hetIndex);
            const double errorLhood =
                    sampleLhood[parentIndex[0]][noiseState] +
                    sampleLhood[parentIndex[1]][noiseState] +
                    sampleLhood[probandIndex][noiseState];

            using namespace TRANSMISSION_STATE;
            stateLhood[ERROR] = log_sum(stateLhood[ERROR], errorLhood);
        }
    }

    std::array<double,TRANSMISSION_STATE::SIZE> statePprob;
    for (unsigned tstate(0); tstate<TRANSMISSION_STATE::SIZE; ++tstate)
    {
        statePprob[tstate] = stateLhood[tstate] + TRANSMISSION_STATE::getPrior(static_cast<TRANSMISSION_STATE::index_t>(tstate));
#ifdef DENOVO_SNV_DEBUG
        const TRANSMISSION_STATE::index_t tidx(static_cast<TRANSMISSION_STATE::index_t>(tstate));
        log_os << "denovo state pprob/lhood/prior: " << TRANSMISSION_STATE::getLabel(tidx)
               << " " << statePprob[tstate] << " " << stateLhood[tstate]
               << " " << TRANSMISSION_STATE::getPrior(tidx) << "\n";
#endif
    }

    //opt_normalize_ln_distro(pprob.begin(),pprob.end(),DDIINDEL_GRID::is_nonsom.val.begin(),rs.max_gt);
    normalize_ln_distro(statePprob.begin(),statePprob.end(),rs.max_gt);

#ifdef DEBUG_INDEL_CALL
    log_os << "INDEL_CALL pprob(noindel),pprob(hom),pprob(het): " << pprob[STAR_DIINDEL::NOINDEL] << " " << pprob[STAR_DIINDEL::HOM] << " " << pprob[STAR_DIINDEL::HET] << "\n";
#endif
    rs.snv_qphred=error_prob_to_qphred(statePprob[TRANSMISSION_STATE::INHERITED] + statePprob[TRANSMISSION_STATE::ERROR]);

}




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

        const auto& sampleCpi(pileups[tierIndex]);

        for (unsigned sampleIndex(0);sampleIndex<sampleSize;++sampleIndex)
        {
            // get likelihoods for each sample
            const CleanedPileup& cpi(*sampleCpi[sampleIndex]);
            const snp_pos_info& pi(cpi.cleanedPileup());
            blt_float_t* lhood(sampleLhood[sampleIndex].data());
            get_diploid_gt_lhood_cached(opt, pi, lhood);
            get_diploid_het_grid_lhood_cached(pi, DIGT_DGRID::HET_RES, lhood);

        }

        calculate_result_set_grid(sinfo, sampleLhood, tier_rs[tierIndex]);
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

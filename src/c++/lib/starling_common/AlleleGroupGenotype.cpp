//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2018 Illumina, Inc.
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

#include "AlleleGroupGenotype.hh"

#include "blt_util/math_util.hh"
#include "blt_util/prob_util.hh"
#include "blt_util/qscore.hh"
#include "blt_util/log.hh"
#include "calibration/IndelErrorModelJson.hh"
#include "common/Exceptions.hh"
#include "starling_common/OrthogonalVariantAlleleCandidateGroupUtil.hh"
#include "starling_common/readMappingAdjustmentUtil.hh"
#include "starling_common/starling_indel_call_pprob_digt.hh"



static
void
updateGenotypeLogLhoodFromAlleleLogLhood(
    const starling_base_deriv_options& dopt,
    const starling_sample_options& sampleOptions,
    const unsigned callerPloidy,
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    const std::vector<double>& alleleLogLhood,
    const ReadPathScores& readScore,
    std::vector<double>& genotypeLogLhood)
{
    static const bool isTier2Pass(false);

    const uint8_t nonRefAlleleCount(alleleGroup.size());
    const uint8_t fullAlleleCount(nonRefAlleleCount+1);

    if (callerPloidy == 1)
    {
        for (unsigned allele0Index(0); allele0Index < fullAlleleCount; ++allele0Index)
        {
            const unsigned genotypeIndex(VcfGenotypeUtil::getGenotypeIndex(allele0Index));
            genotypeLogLhood[genotypeIndex] +=
                integrateOutMappingStatus(dopt, readScore.nonAmbiguousBasesInRead, alleleLogLhood[allele0Index],
                                          isTier2Pass);
        }
    }
    else if (callerPloidy == 2)
    {
        for (unsigned allele1Index(0); allele1Index < fullAlleleCount; ++allele1Index)
        {
            for (unsigned allele0Index(0); allele0Index <= allele1Index; ++allele0Index)
            {
                const bool isHet(allele0Index != allele1Index);
                const unsigned genotypeIndex(VcfGenotypeUtil::getGenotypeIndex(allele0Index, allele1Index));


                double rawLogLhood(0);
                if (isHet)
                {
                    static const double hetAlleleRatio(0.5);
                    static const double loghalf(std::log(0.5));

                    assert(allele1Index > 0);
                    double logHetAllele0Prior(loghalf);
                    double logHetAllele1Prior(loghalf);
                    const IndelKey& allele1Key(alleleGroup.key(allele1Index - 1));
                    get_het_observed_allele_ratio(readScore.read_length, sampleOptions.min_read_bp_flank,
                                                  allele1Key, hetAlleleRatio, logHetAllele0Prior, logHetAllele1Prior);

                    if (allele0Index > 0)
                    {
                        // het-alt genotype -- approximate the expected allele ratio in this case
                        double logRefPrior(loghalf);
                        logHetAllele0Prior = loghalf;
                        const IndelKey& allele0Key(alleleGroup.key(allele0Index - 1));
                        get_het_observed_allele_ratio(readScore.read_length, sampleOptions.min_read_bp_flank,
                                                      allele0Key, hetAlleleRatio, logRefPrior, logHetAllele0Prior);

                        const double normalizeHetRatio(getLogSum(logHetAllele0Prior, logHetAllele1Prior));
                        logHetAllele0Prior -= normalizeHetRatio;
                        logHetAllele1Prior -= normalizeHetRatio;
                    }

                    rawLogLhood = getLogSum(alleleLogLhood[allele0Index] + logHetAllele0Prior,
                                            alleleLogLhood[allele1Index] + logHetAllele1Prior);
                }
                else
                {
                    rawLogLhood = alleleLogLhood[allele0Index];
                }

                genotypeLogLhood[genotypeIndex] +=
                    integrateOutMappingStatus(dopt, readScore.nonAmbiguousBasesInRead, rawLogLhood, isTier2Pass);
            }
        }
    }
    else
    {
        assert(false and "Unexpected ploidy value");
    }
}



///
/// \param[in,out] alleleLhood allele likelihoods, note that these are altered to integrate out mapping status and
///                            be normalized in this func, 0 index is reference
/// \param locusReadStats
///
static
void
updateSupportingReadStats(
    const starling_base_deriv_options& dopt,
    const double readSupportThreshold,
    const uint16_t nonAmbiguousBasesInRead,
    const bool isFwdStrand,
    std::vector<double>& alleleLoglhoods,
    LocusSupportingReadStats& locusReadStats)
{
    static const bool isTier2Pass(false);
    for (auto& alleleHood : alleleLoglhoods)
    {
        alleleHood = integrateOutMappingStatus(dopt, nonAmbiguousBasesInRead, alleleHood, isTier2Pass);
    }
    unsigned maxIndex(0);
    normalizeLogDistro(alleleLoglhoods.begin(), alleleLoglhoods.end(), maxIndex);

    bool isConfidentAlleleFound(false);
    const unsigned fullAlleleCount(alleleLoglhoods.size());
    for (unsigned alleleIndex(0); alleleIndex<fullAlleleCount; ++alleleIndex)
    {
        if (alleleLoglhoods[alleleIndex] < readSupportThreshold) continue;
        locusReadStats.getCounts(isFwdStrand).incrementAlleleCount(alleleIndex);
        isConfidentAlleleFound=true;
        break;
    }

    if (not isConfidentAlleleFound)
    {
        locusReadStats.getCounts(isFwdStrand).nonConfidentCount++;
    }
}


static
const ReadPathScores&
getExemplarReadScore(
    const unsigned sampleIndex,
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    const unsigned readId)
{
    const uint8_t nonRefAlleleCount(alleleGroup.size());

    const ReadPathScores* readScorePtr(nullptr);
    for (unsigned nonRefAlleleIndex(0); nonRefAlleleIndex < nonRefAlleleCount; nonRefAlleleIndex++)
    {
        const IndelSampleData& isd(alleleGroup.data(nonRefAlleleIndex).getSampleData(sampleIndex));
        const auto iditer(isd.read_path_lnp.find(readId));
        if (iditer != isd.read_path_lnp.end())
        {
            readScorePtr = &(iditer->second);
            break;
        }
    }
    assert(nullptr != readScorePtr);
    return (*readScorePtr);
}



void
getVariantAlleleGroupGenotypeLhoodsForSample(
    const starling_base_options& opt,
    const starling_base_deriv_options& dopt,
    const starling_sample_options& sampleOptions,
    const unsigned callerPloidy,
    const unsigned sampleIndex,
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    const OrthogonalVariantAlleleCandidateGroup& contrastGroup,
    std::vector<double>& genotypeLogLhood,
    LocusSupportingReadStats& locusReadStats)
{
    assert(callerPloidy>0u);
    assert(callerPloidy<3u);

    const uint8_t nonRefAlleleCount(alleleGroup.size());
    const uint8_t fullAlleleCount(nonRefAlleleCount+1);
    const unsigned genotypeCount(VcfGenotypeUtil::getGenotypeCount(callerPloidy, fullAlleleCount));

    genotypeLogLhood.resize(genotypeCount);
    std::fill(genotypeLogLhood.begin(), genotypeLogLhood.end(), 0.);

    //---------------------------------------------------
    // get genotype likelihoods
    //
    if (nonRefAlleleCount>0)
    {
        // threshold used to generate supporting count summary (not used for GT likelihoods):
        const double readSupportTheshold(opt.readConfidentSupportThreshold.numval());
        locusReadStats.setAltCount(nonRefAlleleCount);

        static const bool isTier1Only(true);
        std::set<unsigned> readIds;
        getAlleleGroupSupportingReadIds(sampleIndex, alleleGroup, readIds, isTier1Only);

        // contrast group contains alleles intended for an "other" category, such as reported by the <*> allele in
        // vcf
        OrthogonalVariantAlleleCandidateGroup extendedAlleleGroup(alleleGroup);
        for (const auto& contrastAlleleIter : contrastGroup.alleles)
        {
            extendedAlleleGroup.alleles.push_back(contrastAlleleIter);
        }
        const uint8_t extendedNonRefAlleleCount(extendedAlleleGroup.size());
        const uint8_t extendedFullAlleleCount(extendedNonRefAlleleCount+1);

        for (const auto readId : readIds)
        {
            std::vector<double> alleleLogLhood(extendedFullAlleleCount);
            getAlleleLogLhoodFromRead(sampleIndex, extendedAlleleGroup, readId, alleleLogLhood);

            if (fullAlleleCount != extendedFullAlleleCount)
            {
                // TEMPORARY: for now, any contrast allele scores are maxed down into the reference, b/c we don't
                // have a way to report them in the output VCF:
                for (unsigned alleleIndex(fullAlleleCount); alleleIndex<extendedFullAlleleCount; ++alleleIndex)
                {
                    if (alleleLogLhood[alleleIndex] > alleleLogLhood[0])
                    {
                        alleleLogLhood[0] = alleleLogLhood[alleleIndex];
                    }
                }
                alleleLogLhood.resize(fullAlleleCount);
            }

            // get an exemplar read score object, doesn't really matter from which allele...
            const ReadPathScores& readScore(getExemplarReadScore(sampleIndex, alleleGroup, readId));

            updateGenotypeLogLhoodFromAlleleLogLhood(dopt, sampleOptions, callerPloidy, alleleGroup, alleleLogLhood,
                                                     readScore, genotypeLogLhood);

            updateSupportingReadStats(
                dopt, readSupportTheshold, readScore.nonAmbiguousBasesInRead, readScore.is_fwd_strand, alleleLogLhood, locusReadStats);
        }
    }
}

GenotypePriorSet::
GenotypePriorSet(
    const std::string& thetaFilename)
{
    std::map<unsigned, std::vector<double> > thetas;
    if (thetaFilename.empty())
    {
        log_os << "WARNING: theta parameter file was not given. Using internal theta values." << "\n";
        thetas = initializeThetas();
    }
    else
    {
        IndelErrorModelParser::importThetaJsonFile(thetaFilename, thetas);
    }
    initializePriors(thetas);
}

std::map<unsigned, std::vector<double> >
GenotypePriorSet::
initializeThetas()
{
    static const unsigned highHpolRepeatCount(16);
    static const std::vector<double> hpolTheta(
    {
        0.000120268,
        5.97777E-05,
        0.000124648,
        0.000260759,
        0.000589544,
        0.002394583,
        0.007417864,
        0.022660355,
        0.04670561,
        0.082031233,
        0.124548518,
        0.149765438,
        0.168051826,
        0.187346626,
        0.207339703,
        0.225843098,
        0.248849306,
        0.27106361,
        0.334718891,
        0.348811678
    });

    static const unsigned hpolThetaSize = hpolTheta.size();
    assert(hpolThetaSize >= highHpolRepeatCount);

    static const unsigned highDinucRepeatCount(9);
    static const std::vector<double> dinucTheta(
    {
        0.000120268,
        8.73757E-05,
        0.000479319,
        0.002678401,
        0.012194565,
        0.03162284,
        0.060846617,
        0.108263861,
        0.163510548,
        0.204456064,
        0.23462438,
        0.267919304,
        0.290588942,
        0.355588567,
        0.369478351,
        0.378290471,
        0.38555006,
        0.393439865,
        0.395844077,
        0.4
    });

    std::map<unsigned, std::vector<double> > thetas;
    thetas[1] = hpolTheta;
    thetas[2] = dinucTheta;

    static const unsigned dinucThetaSize = dinucTheta.size();
    assert(dinucThetaSize >= highDinucRepeatCount);

    return thetas;
}

void
GenotypePriorSet::
initializePriors(
    const std::map<unsigned, std::vector<double> >& thetas)
{
    static const unsigned maxRepeatingPatternSize(thetas.size());
    assert(maxRepeatingPatternSize == 2);

    _priors.resize(maxRepeatingPatternSize);
    for (unsigned repeatingPatternSize(1); repeatingPatternSize <= maxRepeatingPatternSize; ++repeatingPatternSize)
    {
        const unsigned repeatingPatternSizeIndex(repeatingPatternSize-1);
        auto& strPatternPriors(_priors[repeatingPatternSizeIndex]);

        const unsigned highSTRRepeatCount = thetas.at(repeatingPatternSize).size();

        strPatternPriors.resize(highSTRRepeatCount);
        for (unsigned patternRepeatCount(1); patternRepeatCount <= highSTRRepeatCount; ++patternRepeatCount)
        {
            const unsigned patternRepeatCountIndex(patternRepeatCount - 1);
            const double theta(thetas.at(repeatingPatternSize)[patternRepeatCountIndex]);
            strPatternPriors[patternRepeatCountIndex].initialize(theta);
        }

    }
}

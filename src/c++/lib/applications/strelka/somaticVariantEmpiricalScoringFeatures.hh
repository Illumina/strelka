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

/*
 *      Author: Morten Kallberg
 */

#pragma once

#include "calibration/featuresetUtil.hh"
#include "calibration/VariantScoringModelBase.hh"

#include <cassert>


struct SOMATIC_SNV_SCORING_FEATURES : public FeatureSet
{
    static
    const FeatureSet&
    getInstance()
    {
        static const auto featureSet = SOMATIC_SNV_SCORING_FEATURES();
        return featureSet;
    }

    static
    const char*
    get_name()
    {
        return "SOMATIC_SNV_SCORING_FEATURES";
    }

    /** any change here must be synchronized with:
        src/python/scoringModelTraining/somatic/lib/features/StrelkaSNV.py
     */
    enum index_t
    {
        SomaticSNVQualityAndHomRefGermlineGenotype,
        NormalSampleRelativeTotalLocusDepth,
        TumorSampleAltAlleleFraction,
        RMSMappingQuality,
        ZeroMappingQualityFraction,
        TumorSampleStrandBias,
        TumorSampleReadPosRankSum,
        AlleleCountLogOddsRatio,
        NormalSampleFilteredDepthFraction,
        TumorSampleFilteredDepthFraction,
        SIZE
    };

    static
    const char*
    get_feature_label(const unsigned idx)
    {
        switch (idx)
        {
        case SomaticSNVQualityAndHomRefGermlineGenotype:
            return "SomaticSNVQualityAndHomRefGermlineGenotype";
        case TumorSampleAltAlleleFraction:
            return "TumorSampleAltAlleleFraction";
        case RMSMappingQuality:
            return "RMSMappingQuality";
        case ZeroMappingQualityFraction:
            return "ZeroMappingQualityFraction";
        case TumorSampleStrandBias:
            return "TumorSampleStrandBias";
        case TumorSampleReadPosRankSum:
            return "TumorSampleReadPosRankSum";
        case NormalSampleRelativeTotalLocusDepth:
            return "NormalSampleRelativeTotalLocusDepth";
        case AlleleCountLogOddsRatio:
            return "AlleleCountLogOddsRatio";
        case NormalSampleFilteredDepthFraction:
            return "NormalSampleFilteredDepthFraction";
        case TumorSampleFilteredDepthFraction:
            return "TumorSampleFilteredDepthFraction";
        default:
            assert(false && "Unknown feature");
            return nullptr;
        }
    }

    const char*
    getName() const override
    {
        return get_name();
    }

    unsigned
    size() const override
    {
        return SIZE;
    }


    const char*
    getFeatureLabel(const unsigned idx) const override
    {
        return get_feature_label(idx);
    }

};


///
/// additional experimental features not used in the current scoring model
///
/// these should only be output as part of a non-default training mode
///
struct SOMATIC_SNV_SCORING_DEVELOPMENT_FEATURES : public FeatureSet
{
    static
    const FeatureSet&
    getInstance()
    {
        static const auto featureSet = SOMATIC_SNV_SCORING_DEVELOPMENT_FEATURES();
        return featureSet;
    }

    static
    const char*
    get_name()
    {
        return "SOMATIC_SNV_SCORING_DEVELOPMENT_FEATURES";
    }

    enum index_t
    {
        TumorSampleAltAlleleMedianReadPosVariation,
        TumorSampleAltAlleleMedianReadPos,
        NormalSampleSpanningDeletionFraction,
        TumorSampleSpanningDeletionFraction,
        SIZE
    };

    static
    const char*
    get_feature_label(const unsigned idx)
    {
        switch (idx)
        {
        case TumorSampleAltAlleleMedianReadPosVariation:
            return "TumorSampleAltAlleleMedianReadPosVariation";
        case TumorSampleAltAlleleMedianReadPos:
            return "TumorSampleAltAlleleMedianReadPos";
        case NormalSampleSpanningDeletionFraction:
            return "NormalSampleSpanningDeletionFraction";
        case TumorSampleSpanningDeletionFraction:
            return "TumorSampleSpanningDeletionFraction";
        default:
            assert(false && "Unknown feature");
            return nullptr;
        }
    }

    const char*
    getName() const override
    {
        return get_name();
    }

    unsigned
    size() const override
    {
        return SIZE;
    }


    const char*
    getFeatureLabel(const unsigned idx) const override
    {
        return get_feature_label(idx);
    }

};



struct SOMATIC_INDEL_SCORING_FEATURES : public FeatureSet
{
    static
    const FeatureSet&
    getInstance()
    {
        static const auto featureSet = SOMATIC_INDEL_SCORING_FEATURES();
        return featureSet;
    }

    const char*
    getName() const override
    {
        return "SOMATIC_INDEL_SCORING_FEATURES";
    }

    /** Make sure the features are the same as used in the model
     */
    enum index_t
    {
        SomaticIndelQualityAndHomRefGermlineGenotype,
        TumorSampleReadPosRankSum,
        TumorSampleLogSymmetricStrandOddsRatio,
        IndelRepeatCount,
        InterruptedHomopolymerLength,
        RefRepeatCount,
        RepeatUnitLength,
        TumorSampleIndelNoiseLogOdds,
        TumorNormalIndelAlleleLogOdds,
        AlleleCountLogOddsRatio,
        SIZE
    };

    unsigned
    size() const override
    {
        return SIZE;
    }


    const char*
    getFeatureLabel(const unsigned idx) const override
    {
        switch (idx)
        {
        case SomaticIndelQualityAndHomRefGermlineGenotype:
            return "SomaticIndelQualityAndHomRefGermlineGenotype";
        case TumorSampleReadPosRankSum:
            return "TumorSampleReadPosRankSum";
        case TumorSampleLogSymmetricStrandOddsRatio:
            return "TumorSampleLogSymmetricStrandOddsRatio";
        case IndelRepeatCount:
            return "IndelRepeatCount";
        case InterruptedHomopolymerLength:
            return "InterruptedHomopolymerLength";
        case RefRepeatCount:
            return "RefRepeatCount";
        case RepeatUnitLength:
            return "RepeatUnitLength";
        case TumorSampleIndelNoiseLogOdds:
            return "TumorSampleIndelNoiseLogOdds";
        case TumorNormalIndelAlleleLogOdds:
            return "TumorNormalIndelAlleleLogOdds";
        case AlleleCountLogOddsRatio:
            return "AlleleCountLogOddsRatio";
        default:
            assert(false && "Unknown feature");
            return nullptr;
        }
    }

};


///
/// additional experimental features not used in the current scoring model
///
/// these should only be output as part of a non-default training mode
///
struct SOMATIC_INDEL_SCORING_DEVELOPMENT_FEATURES : public FeatureSet
{
    static
    const FeatureSet&
    getInstance()
    {
        static const auto featureSet = SOMATIC_INDEL_SCORING_DEVELOPMENT_FEATURES();
        return featureSet;
    }

    static
    const char*
    get_name()
    {
        return "SOMATIC_INDEL_SCORING_DEVELOPMENT_FEATURES";
    }

    enum index_t
    {
        N_AF,
        T_AF,
        N_OF,
        T_OF,
        N_BCN,
        T_BCN,
        TumorSampleAbsReadPosRankSum,
        TumorSampleLogStrandOddsRatio,
        TumorSampleAbsLogStrandOddsRatio,
        SIZE
    };

    static
    const char*
    get_feature_label(const unsigned idx)
    {
        switch (idx)
        {
        case N_AF:
            return "N_AF";
        case T_AF:
            return "T_AF";
        case N_OF:
            return "N_OF";
        case T_OF:
            return "T_OF";
        case N_BCN:
            return "N_BCN";
        case T_BCN:
            return "T_BCN";
        case TumorSampleAbsReadPosRankSum:
            return "TumorSampleAbsReadPosRankSum";
        case TumorSampleLogStrandOddsRatio:
            return "TumorSampleLogStrandOddsRatio";
        case TumorSampleAbsLogStrandOddsRatio:
            return "TumorSampleAbsLogStrandOddsRatio";
        default:
            assert(false && "Unknown feature");
            return nullptr;
        }
    }

    const char*
    getName() const override
    {
        return get_name();
    }

    unsigned
    size() const override
    {
        return SIZE;
    }


    const char*
    getFeatureLabel(const unsigned idx) const override
    {
        return get_feature_label(idx);
    }

};


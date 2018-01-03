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

///  \author mkallberg
///

#pragma once

#include "calibration/featuresetUtil.hh"

#include <cassert>


struct GERMLINE_SNV_SCORING_FEATURES : public FeatureSet
{
    static
    const FeatureSet&
    getInstance()
    {
        static const auto featureSet = GERMLINE_SNV_SCORING_FEATURES();
        return featureSet;
    }

    const char*
    getName() const override
    {
        return "GERMLINE_SNV_SCORING_FEATURES";
    }

    /** any change here must be done together with changing
        src/python/scoringModelTraining/germline/lib/evs/features/GermlineSNV.py
     */
    enum index_t
    {
        GenotypeCategory,
        SampleRMSMappingQuality,
        SiteHomopolymerLength,
        SampleStrandBias,
        SampleRMSMappingQualityRankSum,
        SampleReadPosRankSum,
        RelativeTotalLocusDepth,
        SampleUsedDepthFraction,
        ConservativeGenotypeQuality,
        NormalizedAltHaplotypeCountRatio,
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
        case GenotypeCategory:
            return "GenotypeCategory";
        case SampleRMSMappingQuality:
            return "SampleRMSMappingQuality";
        case SiteHomopolymerLength:
            return "SiteHomopolymerLength";
        case SampleStrandBias:
            return "SampleStrandBias";
        case SampleRMSMappingQualityRankSum:
            return "SampleRMSMappingQualityRankSum";
        case SampleReadPosRankSum:
            return "SampleReadPosRankSum";
        case RelativeTotalLocusDepth:
            return "RelativeTotalLocusDepth";
        case SampleUsedDepthFraction:
            return "SampleUsedDepthFraction";
        case ConservativeGenotypeQuality:
            return "ConservativeGenotypeQuality";
        case NormalizedAltHaplotypeCountRatio:
            return "NormalizedAltHaplotypeCountRatio";
        default:
            assert(false && "Unknown feature");
            return nullptr;
        }
    }
};


/// additional experimental features not used in the current scoring model
///
/// these should only be output as part of a non-default training mode
///
struct GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES : public FeatureSet
{
    static
    const FeatureSet&
    getInstance()
    {
        static const auto featureSet = GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES();
        return featureSet;
    }

    const char*
    getName() const override
    {
        return "GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES";
    }

    enum index_t
    {
        I_BaseQRankSum,
        SampleIndelAlleleBiasLower,
        SampleIndelAlleleBias,
        I_RawPos,
        I_RawBaseQ,
        mapqZeroFraction,
        QUAL_NORM,
        F_GQX_NORM,
        F_GQ_NORM,
        SampleRefAlleleDepthFraction,
        VariantAlleleQuality,
        F_GQ,
        SamplePrimaryAltAlleleDepthFraction,
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
        case I_BaseQRankSum:
            return "I_BaseQRankSum";
        case SampleIndelAlleleBiasLower:
            return "SampleIndelAlleleBiasLower";
        case SampleIndelAlleleBias:
            return "SampleIndelAlleleBias";
        case I_RawPos:
            return "I_RawPos";
        case I_RawBaseQ:
            return "I_RawBaseQ";
        case mapqZeroFraction:
            return "mapqZeroFraction";
        case QUAL_NORM:
            return "QUAL_NORM";
        case F_GQX_NORM:
            return "F_GQX_NORM";
        case F_GQ_NORM:
            return "F_GQ_NORM";
        case SampleRefAlleleDepthFraction:
            return "SampleRefAlleleDepthFraction";
        case VariantAlleleQuality:
            return "VariantAlleleQuality";
        case F_GQ:
            return "F_GQ";
        case SamplePrimaryAltAlleleDepthFraction:
            return "SamplePrimaryAltAlleleDepthFraction";
        default:
            assert(false && "Unknown feature");
            return nullptr;
        }
    }
};



struct GERMLINE_INDEL_SCORING_FEATURES : public FeatureSet
{
    static
    const FeatureSet&
    getInstance()
    {
        static const auto featureSet = GERMLINE_INDEL_SCORING_FEATURES();
        return featureSet;
    }

    const char*
    getName() const override
    {
        return "GERMLINE_INDEL_SCORING_FEATURES";
    }

    /** Make sure the features are the same as used in the model.
        Any change here must be done together with changing
        src/python/scoringModelTraining/germline/lib/evs/features/GermlineIndel.py
    */
    enum index_t
    {
        GenotypeCategory,
        SampleIndelRepeatCount,
        SampleIndelRepeatUnitSize,
        SampleIndelAlleleBiasLower,
        SampleIndelAlleleBias,
        SampleProxyRMSMappingQuality,
        RelativeTotalLocusDepth,
        SamplePrimaryAltAlleleDepthFraction,
        ConservativeGenotypeQuality,
        InterruptedHomopolymerLength,
        ContextCompressability,
        IndelCategory,
        NormalizedAltHaplotypeCountRatio,
        SampleAlleleCountStrandBias,
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
        case GenotypeCategory:
            return "GenotypeCategory";
        case SampleIndelRepeatCount:
            return "SampleIndelRepeatCount";
        case SampleIndelRepeatUnitSize:
            return "SampleIndelRepeatUnitSize";
        case SampleIndelAlleleBiasLower:
            return "SampleIndelAlleleBiasLower";
        case SampleIndelAlleleBias:
            return "SampleIndelAlleleBias";
        case SampleProxyRMSMappingQuality:
            return "SampleProxyRMSMappingQuality";
        case RelativeTotalLocusDepth:
            return "RelativeTotalLocusDepth";
        case SamplePrimaryAltAlleleDepthFraction:
            return "SamplePrimaryAltAlleleDepthFraction";
        case ConservativeGenotypeQuality:
            return "ConservativeGenotypeQuality";
        case InterruptedHomopolymerLength:
            return "InterruptedHomopolymerLength";
        case ContextCompressability:
            return "ContextCompressability";
        case IndelCategory:
            return "IndelCategory";
        case NormalizedAltHaplotypeCountRatio:
            return "NormalizedAltHaplotypeCountRatio";
        case SampleAlleleCountStrandBias:
            return "SampleAlleleCountStrandBias";
        default:
            assert(false && "Unknown feature");
            return nullptr;
        }
    }
};


/// additional experimental features not used in the current scoring model
///
/// these should only be output as part of a non-default training mode
///
struct GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES : public FeatureSet
{
    static
    const FeatureSet&
    getInstance()
    {
        static const auto featureSet = GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES();
        return featureSet;
    }

    const char*
    getName() const override
    {
        return "GERMLINE_INDEL_DEVELOPMENT_SCORING_FEATURES";
    }

    enum index_t
    {
        SampleRefRepeatCount,
        mapqZeroFraction,
        F_DPI_NORM,
        QUAL_NORM,
        F_GQX_NORM,
        F_GQ_NORM,
        SampleRefAlleleDepthFraction,
        VariantAlleleQuality,
        F_GQ,
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
        case SampleRefRepeatCount:
            return "SampleRefRepeatCount";
        case mapqZeroFraction:
            return "mapqZeroFraction";
        case F_DPI_NORM:
            return "F_DPI_NORM";
        case QUAL_NORM:
            return "QUAL_NORM";
        case F_GQX_NORM:
            return "F_GQX_NORM";
        case F_GQ_NORM:
            return "F_GQ_NORM";
        case SampleRefAlleleDepthFraction:
            return "SampleRefAlleleDepthFraction";
        case VariantAlleleQuality:
            return "VariantAlleleQuality";
        case F_GQ:
            return "F_GQ";
        default:
            assert(false && "Unknown feature");
            return nullptr;
        }
    }
};

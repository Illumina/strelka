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
 *  \author mkallberg
 */

#pragma once

#include "calibration/featuresetUtil.hh"

#include <cassert>


struct RNA_SNV_SCORING_FEATURES : public FeatureSet
{
    static
    const FeatureSet&
    getInstance()
    {
        static const auto featureSet = RNA_SNV_SCORING_FEATURES();
        return featureSet;
    }

    const char*
    getName() const override
    {
        return "RNA_SNV_SCORING_FEATURES";
    }

    /** any change here must be done together with changing
        src/python/scoringModelTraining/rna/lib/features/StrelkaSNV.py
     */
    enum index_t
    {

        SiteHomopolymerLength,
        SampleStrandBias,
        SamplePrimaryAltAlleleDepth,
        VariantAlleleQuality,
        SampleMeanDistanceFromReadEdge,
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
        case SiteHomopolymerLength:
            return "SiteHomopolymerLength";
        case SampleStrandBias:
            return "SampleStrandBias";
        case SamplePrimaryAltAlleleDepth:
            return "SamplePrimaryAltAlleleDepth";
        case SamplePrimaryAltAlleleDepthFraction:
            return "SamplePrimaryAltAlleleDepthFraction";
        case VariantAlleleQuality:
            return "VariantAlleleQuality";
        case SampleMeanDistanceFromReadEdge:
            return "SampleMeanDistanceFromReadEdge";
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
struct RNA_SNV_SCORING_DEVELOPMENT_FEATURES : public FeatureSet
{
    static
    const FeatureSet&
    getInstance()
    {
        static const auto featureSet = RNA_SNV_SCORING_DEVELOPMENT_FEATURES();
        return featureSet;
    }

    const char*
    getName() const override
    {
        return "RNA_SNV_SCORING_DEVELOPMENT_FEATURES";
    }

    enum index_t
    {
        GT,
        I_AvgBaseQ,
        I_BaseQRankSum,
        F_DP,
        F_DPF,
        SampleReadPosRankSum,
        SampleRefAlleleDepth,
        SampleRMSMappingQuality,
        SampleRMSMappingQualityRankSum,
        mapqZeroFraction,
        SampleUsedDepthFraction,
        QUAL_NORM,
        F_GQX_NORM,
        F_GQ_NORM,
        SampleRefAlleleDepthFraction,
        ConservativeGenotypeQuality,
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
        case GT:
            return "GT";
        case F_DP:
            return "F_DP";
        case F_DPF:
            return "F_DPF";
        case I_AvgBaseQ:
            return "I_AvgBaseQ";
        case I_BaseQRankSum:
            return "I_BaseQRankSum";
        case SampleReadPosRankSum:
            return "SampleReadPosRankSum";
        case SampleRefAlleleDepth:
            return "SampleRefAlleleDepth";
        case SampleRMSMappingQuality:
            return "SampleRMSMappingQuality";
        case SampleRMSMappingQualityRankSum:
            return "SampleRMSMappingQualityRankSum";
        case mapqZeroFraction:
            return "mapqZeroFraction";
        case SampleUsedDepthFraction:
            return "SampleUsedDepthFraction";
        case QUAL_NORM:
            return "QUAL_NORM";
        case F_GQX_NORM:
            return "F_GQX_NORM";
        case F_GQ_NORM:
            return "F_GQ_NORM";
        case SampleRefAlleleDepthFraction:
            return "SampleRefAlleleDepthFraction";
        case ConservativeGenotypeQuality:
            return "ConservativeGenotypeQuality";
        case F_GQ:
            return "F_GQ";
        default:
            assert(false && "Unknown feature");
            return nullptr;
        }
    }
};



struct RNA_INDEL_SCORING_FEATURES : public FeatureSet
{
    static
    const FeatureSet&
    getInstance()
    {
        static const auto featureSet = RNA_INDEL_SCORING_FEATURES();
        return featureSet;
    }

    const char*
    getName() const override
    {
        return "RNA_INDEL_SCORING_FEATURES";
    }

    /** Make sure the features are the same as used in the python training code
     */
    enum index_t
    {
        SampleRefAlleleDepth,
        SamplePrimaryAltAlleleDepth,
        SampleIndelRepeatCount,
        SampleIndelRepeatUnitSize,
        VariantAlleleQuality,
        SampleIndelMeanDistanceFromReadEdge,
        SampleRefRepeatCount,
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
        case SampleRefRepeatCount:
            return "SampleRefRepeatCount";
        case SampleIndelRepeatCount:
            return "SampleIndelRepeatCount";
        case SampleIndelRepeatUnitSize:
            return "SampleIndelRepeatUnitSize";
        case SampleRefAlleleDepth:
            return "SampleRefAlleleDepth";
        case SamplePrimaryAltAlleleDepth:
            return "SamplePrimaryAltAlleleDepth";
        case SampleIndelMeanDistanceFromReadEdge:
            return "SampleIndelMeanDistanceFromReadEdge";
        case SamplePrimaryAltAlleleDepthFraction:
            return "SamplePrimaryAltAlleleDepthFraction";
        case VariantAlleleQuality:
            return "VariantAlleleQuality";
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
struct RNA_INDEL_SCORING_DEVELOPMENT_FEATURES : public FeatureSet
{
    static
    const FeatureSet&
    getInstance()
    {
        static const auto featureSet = RNA_INDEL_SCORING_DEVELOPMENT_FEATURES();
        return featureSet;
    }

    const char*
    getName() const override
    {
        return "RNA_INDEL_DEVELOPMENT_SCORING_FEATURES";
    }

    enum index_t
    {
        SampleIndelAlleleBiasLower,
        SampleIndelAlleleBias,
        F_DPI,
        SampleProxyRMSMappingQuality,
        mapqZeroFraction,
        F_DPI_NORM,
        QUAL_NORM,
        F_GQX_NORM,
        F_GQ_NORM,
        SampleRefAlleleDepthFraction,
        ConservativeGenotypeQuality,
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
        case SampleIndelAlleleBiasLower:
            return "SampleIndelAlleleBiasLower";
        case SampleIndelAlleleBias:
            return "SampleIndelAlleleBias";
        case F_DPI:
            return "F_DPI";
        case SampleProxyRMSMappingQuality:
            return "SampleProxyRMSMappingQuality";
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
        case ConservativeGenotypeQuality:
            return "ConservativeGenotypeQuality";
        case F_GQ:
            return "F_GQ";
        default:
            assert(false && "Unknown feature");
            return nullptr;
        }
    }
};

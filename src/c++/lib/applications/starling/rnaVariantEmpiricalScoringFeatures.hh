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
        GT,
        QUAL,
        F_DP,
        F_DPF,
        F_GQ,
        F_GQX,
        I_AvgBaseQ,
        I_AvgPos,
        I_BaseQRankSum,
        SampleReadPosRankSum,
        SiteHomopolymerLength,
        SampleStrandBias,
        AD0,
        AD1,
        ADR,
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
        case QUAL:
            return "QUAL";
        case F_DP:
            return "F_DP";
        case F_DPF:
            return "F_DPF";
        case F_GQ:
            return "F_GQ";
        case F_GQX:
            return "F_GQX";
        case I_AvgBaseQ:
            return "I_AvgBaseQ";
        case I_AvgPos:
            return "I_AvgPos";
        case I_BaseQRankSum:
            return "I_BaseQRankSum";
        case SampleReadPosRankSum:
            return "SampleReadPosRankSum";
        case SiteHomopolymerLength:
            return "SiteHomopolymerLength";
        case SampleStrandBias:
            return "SampleStrandBias";
        case AD0:
            return "AD0";
        case AD1:
            return "AD1";
        case ADR:
            return "ADR";
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
        SampleRMSMappingQuality,
        SampleRMSMappingQualityRankSum,
        mapqZeroFraction,
        SampleUsedDepthFraction,
        QUAL_NORM,
        F_GQX_NORM,
        F_GQ_NORM,
        AD0_NORM,
        SamplePrimaryAltAlleleDepthFraction,
        QUAL_EXACT,
        ConservativeGenotypeQuality,
        F_GQ_EXACT,
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
        case AD0_NORM:
            return "AD0_NORM";
        case SamplePrimaryAltAlleleDepthFraction:
            return "SamplePrimaryAltAlleleDepthFraction";
        case QUAL_EXACT:
            return "QUAL_EXACT";
        case ConservativeGenotypeQuality:
            return "ConservativeGenotypeQuality";
        case F_GQ_EXACT:
            return "F_GQ_EXACT";
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
        QUAL,
        F_GQX,
        REFREP1,
        IndelAlleleRepeatCount,
        IndelAlleleRepeatUnitSize,
        AD0,
        AD1,
        SampleIndelAlleleBiasLower,
        SampleIndelAlleleBias,
        F_DPI,
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
        case QUAL:
            return "QUAL";
        case F_GQX:
            return "F_GQX";
        case REFREP1:
            return "REFREP1";
        case IndelAlleleRepeatCount:
            return "IndelAlleleRepeatCount";
        case IndelAlleleRepeatUnitSize:
            return "IndelAlleleRepeatUnitSize";
        case AD0:
            return "AD0";
        case AD1:
            return "AD1";
        case SampleIndelAlleleBiasLower:
            return "SampleIndelAlleleBiasLower";
        case SampleIndelAlleleBias:
            return "SampleIndelAlleleBias";
        case F_DPI:
            return "F_DPI";
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
        F_GQ,
        SampleProxyRMSMappingQuality,
        mapqZeroFraction,
        F_DPI_NORM,
        QUAL_NORM,
        F_GQX_NORM,
        F_GQ_NORM,
        AD0_NORM,
        SamplePrimaryAltAlleleDepthFraction,
        QUAL_EXACT,
        ConservativeGenotypeQuality,
        F_GQ_EXACT,
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
        case F_GQ:
            return "F_GQ";
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
        case AD0_NORM:
            return "AD0_NORM";
        case SamplePrimaryAltAlleleDepthFraction:
            return "SamplePrimaryAltAlleleDepthFraction";
        case QUAL_EXACT:
            return "QUAL_EXACT";
        case ConservativeGenotypeQuality:
            return "ConservativeGenotypeQuality";
        case F_GQ_EXACT:
            return "F_GQ_EXACT";
        default:
            assert(false && "Unknown feature");
            return nullptr;
        }
    }
};

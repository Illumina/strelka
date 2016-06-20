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
        static const RNA_SNV_SCORING_FEATURES featureSet;
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
        I_ReadPosRankSum,
        I_SNVHPOL,
        I_SNVSB,
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
        case I_ReadPosRankSum:
            return "I_ReadPosRankSum";
        case I_SNVHPOL:
            return "I_SNVHPOL";
        case I_SNVSB:
            return "I_SNVSB";
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
        static const RNA_SNV_SCORING_DEVELOPMENT_FEATURES featureSet;
        return featureSet;
    }

    const char*
    getName() const override
    {
        return "RNA_SNV_SCORING_DEVELOPMENT_FEATURES";
    }

    enum index_t
    {
        I_MQ,
        I_MQRankSum,
        mapqZeroFraction,
        F_DP_NORM,
        QUAL_NORM,
        F_GQX_NORM,
        F_GQ_NORM,
        AD0_NORM,
        AD1_NORM,
        QUAL_EXACT,
        F_GQX_EXACT,
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
        case I_MQ:
            return "I_MQ";
        case I_MQRankSum:
            return "I_MQRankSum";
        case mapqZeroFraction:
            return "mapqZeroFraction";
        case F_DP_NORM:
            return "F_DP_NORM";
        case QUAL_NORM:
            return "QUAL_NORM";
        case F_GQX_NORM:
            return "F_GQX_NORM";
        case F_GQ_NORM:
            return "F_GQ_NORM";
        case AD0_NORM:
            return "AD0_NORM";
        case AD1_NORM:
            return "AD1_NORM";
        case QUAL_EXACT:
            return "QUAL_EXACT";
        case F_GQX_EXACT:
            return "F_GQX_EXACT";
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
        static const RNA_INDEL_SCORING_FEATURES featureSet;
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
        IDREP1,
        RULEN1,
        AD0,
        AD1,
        AD2,
        ABlower,
        AB,
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
        case IDREP1:
            return "IDREP1";
        case RULEN1:
            return "RULEN1";
        case AD0:
            return "AD0";
        case AD1:
            return "AD1";
        case AD2:
            return "AD2";
        case ABlower:
            return "ABlower";
        case AB:
            return "AB";
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
        static const RNA_INDEL_SCORING_DEVELOPMENT_FEATURES featureSet;
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
        F_MQ,
        mapqZeroFraction,
        F_DPI_NORM,
        QUAL_NORM,
        F_GQX_NORM,
        F_GQ_NORM,
        AD0_NORM,
        AD1_NORM,
        AD2_NORM,
        QUAL_EXACT,
        F_GQX_EXACT,
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
        case F_MQ:
            return "F_MQ";
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
        case AD1_NORM:
            return "AD1_NORM";
        case AD2_NORM:
            return "AD2_NORM";
        case QUAL_EXACT:
            return "QUAL_EXACT";
        case F_GQX_EXACT:
            return "F_GQX_EXACT";
        case F_GQ_EXACT:
            return "F_GQ_EXACT";
        default:
            assert(false && "Unknown feature");
            return nullptr;
        }
    }
};

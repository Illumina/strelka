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
 *      Author: mkallberg
 */

#pragma once

#include "calibration/featuresetUtil.hh"
#include "calibration/VariantScoringModelBase.hh"

#include <cassert>


struct GERMLINE_SNV_SCORING_FEATURES
{
    static
    const char*
    get_name()
    {
        return "GERMLINE_SNV_SCORING_FEATURES";
    }

    /** any change here must be done together with changing
        src/python/scoringModelTraining/germline/lib/features/StrelkaSNV.py
     */
    enum index_t
    {
        QUAL,
        F_GQX,
        F_GQ,
        I_SNVSB,
        I_SNVHPOL,
        F_DP,
        F_DPF,
        AD0,
        AD1,
        I_MQ,
        I_ReadPosRankSum,
        I_BaseQRankSum,
        I_MQRankSum,
        ABlower,
        AB,
        SIZE
    };

    static
    const char*
    get_feature_label(const unsigned idx)
    {
        switch (idx)
        {
        case QUAL:
            return "QUAL";
        case F_GQX:
            return "F_GQX";
        case F_GQ:
            return "F_GQ";
        case I_SNVSB:
            return "I_SNVSB";
        case I_SNVHPOL:
            return "I_SNVHPOL";
        case F_DP:
            return "F_DP";
        case F_DPF:
            return "F_DPF";
        case AD0:
            return "AD0";
        case AD1:
            return "AD1";
        case I_MQ:
            return "I_MQ";
        case I_ReadPosRankSum:
            return "I_ReadPosRankSum";
        case I_BaseQRankSum:
            return "I_BaseQRankSum";
        case I_MQRankSum:
            return "I_MQRankSum";
        case ABlower:
            return "ABlower";
        case AB:
            return "AB";
        default:
            assert(false && "Unknown feature");
            return nullptr;
        }
    }

    static
    const VariantScoringModelMetadata::featureMap_t&
    getFeatureMap()
    {
        static const FeatureMapMaker<GERMLINE_SNV_SCORING_FEATURES> fmm;
        return fmm.result;
    }
};


///
/// additional experimental features not used in the current scoring model
///
/// these should only be output as part of a non-default training mode
///
struct GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES
{
    static
    const char*
    get_name()
    {
        return "GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES";
    }

    enum index_t
    {
        I_RawPos,
        I_RawBaseQ,
        mapqZeroFraction,
        F_DP_NORM,
        TDP_NORM,
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

    static
    const char*
    get_feature_label(const unsigned idx)
    {
        switch (idx)
        {
        case I_RawPos:
            return "I_RawPos";
        case I_RawBaseQ:
            return "I_RawBaseQ";
        case mapqZeroFraction:
            return "mapqZeroFraction";
        case F_DP_NORM:
            return "F_DP_NORM";
        case TDP_NORM:
            return "TDP_NORM";
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



struct GERMLINE_INDEL_SCORING_FEATURES
{
    static
    const char*
    get_name()
    {
        return "GERMLINE_INDEL_SCORING_FEATURES";
    }

    /** Make sure the features are the same as used in the model
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

    static
    const char*
    get_feature_label(const unsigned idx)
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

    static
    const VariantScoringModelMetadata::featureMap_t&
    getFeatureMap()
    {
        static const FeatureMapMaker<GERMLINE_INDEL_SCORING_FEATURES> fmm;
        return fmm.result;
    }
};


/// additional experimental features not used in the current scoring model
///
/// these should only be output as part of a non-default training mode
///
struct GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES
{
    static
    const char*
    get_name()
    {
        return "GERMLINE_INDEL_DEVELOPMENT_SCORING_FEATURES";
    }

    enum index_t
    {
        F_GQ,
        F_MQ,
        mapqZeroFraction,
        F_DPI_NORM,
        TDP_NORM,
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

    static
    const char*
    get_feature_label(const unsigned idx)
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
        case TDP_NORM:
            return "TDP_NORM";
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

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
    enum index_t
    {
        I_RawPos,
        I_RawBaseQ,
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
        default:
            assert(false && "Unknown feature");
            return nullptr;
        }
    }
};



struct GERMLINE_INDEL_SCORING_FEATURES
{
    /** Make sure the features are the same as used in the model
     */
    enum index_t
    {
        QSI_NT,
        ABS_T_RR,
        ABS_T_SOR,
        IC,
        IHP,
        RC,
        RU_LEN,
        TNR,
        AFR,
        LOR,
        SIZE
    };

    static
    const char*
    get_feature_label(const unsigned idx)
    {
        switch (idx)
        {
        case QSI_NT:
            return "QSI_NT";
        case ABS_T_RR:
            return "ABS_T_RR";
        case ABS_T_SOR:
            return "ABS_T_SOR";
        case IC:
            return "IC";
        case IHP:
            return "IHP";
        case RC:
            return "RC";
        case RU_LEN:
            return "RU_LEN";
        case TNR:
            return "TNR";
        case AFR:
            return "AFR";
        case LOR:
            return "LOR";
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


///
/// additional experimental features not used in the current scoring model
///
/// these should only be output as part of a non-default training mode
///
struct GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES
{
    enum index_t
    {
        N_AF,
        T_AF,
        N_OF,
        T_OF,
        N_BCN,
        T_BCN,
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
        default:
            assert(false && "Unknown feature");
            return nullptr;
        }
    }
};


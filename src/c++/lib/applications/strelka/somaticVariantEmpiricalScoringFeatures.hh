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
        static const SOMATIC_SNV_SCORING_FEATURES featureSet;
        return featureSet;
    }

    static
    const char*
    get_name()
    {
        return "SOMATIC_SNV_SCORING_FEATURES";
    }

    /** any change here must be done together with changing
        src/python/scoringModelTraining/somatic/lib/features/StrelkaSNV.py
     */
    enum index_t
    {
        QSS_NT,
        N_FDP_RATE,
        T_FDP_RATE,
        N_SDP_RATE,
        T_SDP_RATE,
        N_DP_RATE,
        TIER1_ALT_RATE,
        MQ,
        n_mapq0,
        strandBias,
        ReadPosRankSum,
        MQ0_FRAC,
        LOR,
        SIZE
    };

    static
    const char*
    get_feature_label(const unsigned idx)
    {
        switch (idx)
        {
        case QSS_NT:
            return "QSS_NT";
        case N_FDP_RATE:
            return "N_FDP_RATE";
        case T_FDP_RATE:
            return "T_FDP_RATE";
        case N_SDP_RATE:
            return "N_SDP_RATE";
        case T_SDP_RATE:
            return "T_SDP_RATE";
        case N_DP_RATE:
            return "N_DP_RATE";
        case TIER1_ALT_RATE:
            return "TIER1_ALT_RATE";
        case MQ:
            return "MQ";
        case n_mapq0:
            return "n_mapq0";
        case strandBias:
            return "strandBias";
        case ReadPosRankSum:
            return "ReadPosRankSum";
        case MQ0_FRAC:
            return "MQ0_FRAC";
        case LOR:
            return "LogOddsRatio";
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
        static const SOMATIC_SNV_SCORING_DEVELOPMENT_FEATURES featureSet;
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
        altmap,
        altpos,
        SIZE
    };

    static
    const char*
    get_feature_label(const unsigned idx)
    {
        switch (idx)
        {
        case altmap:
            return "altmap";
        case altpos:
            return "altpos";
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
        static const SOMATIC_INDEL_SCORING_FEATURES featureSet;
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
        static const SOMATIC_INDEL_SCORING_DEVELOPMENT_FEATURES featureSet;
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


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

#include <cassert>

struct STRELKA_SNV_VQSR_FEATURES
{

/** any change here must be done together with changing
    src/python/somaticVQSRTraining/lib/features//StrelkaSNV.py
 */
    enum index_t
    {
        QSS_NT,
        N_FDP_RATE,
        T_FDP_RATE,
        N_SDP_RATE,
        T_SDP_RATE,
        N_DP_RATE,
        TIER1_ALLELE_RATE,
        MQ,
        n_mapq0,
        strandBias,
        ReadPosRankSum,
        altmap,
        altpos,
        pnoise,
        pnoise2,
        SIZE
    };

    inline
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
            case TIER1_ALLELE_RATE:
                return "TIER1_ALLELE_RATE";
            case MQ:
                return "MQ";
            case n_mapq0:
                return "n_mapq0";
            case strandBias:
                return "strandBias";
            case ReadPosRankSum:
                return "ReadPosRankSum";
            case altmap:
                return "altmap";
            case altpos:
                return "altpos";
            case pnoise:
                return "pnoise";
            case pnoise2:
                return "pnoise2";
            default:
                assert(false && "Unknown feature");
                return nullptr;
        }
    }
};

struct STRELKA_INDEL_VQSR_FEATURES
{

/** Make sure the features are the same as used in the model
 */
    enum index_t
    {
        QSI_NT,
        AOR,
        OD,
        ABS_T_RR,
        T_FS,
        ABS_T_SOR,
        SIZE
    };

    inline
    const char*
    get_feature_label(const unsigned idx)
    {
        switch (idx)
        {
            case QSI_NT:
                return "QSI_NT";
            case AOR:
                return "AOR";
            case OD:
                return "OD";
            case ABS_T_RR:
                return "ABS_T_RR";
            case T_FS:
                return "T_FS";
            case ABS_T_SOR:
                return "ABS_T_SOR";
            default:
                assert(false && "Unknown feature");
                return nullptr;
        }
    }
};

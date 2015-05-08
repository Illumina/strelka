// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Starka
// Copyright (c) 2009-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//
/*
 *      Author: mkallberg
 */

#pragma once

#include <cassert>

namespace STRELKA_VQSR_FEATURES
{

/** any change here must be done together with changing
    src/python/lib/vqsr/StrelkaSNV.py
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
}

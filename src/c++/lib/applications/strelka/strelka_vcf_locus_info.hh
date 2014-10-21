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

///
/// \author Chris Saunders
///

#pragma once

#include <cassert>

#include <bitset>
#include <iosfwd>


namespace STRELKA_VCF_FILTERS
{

enum index_t
{
    // SNVs and indels:
    HighDepth,
    // SNVs only:
    BCNoise,
    SpanDel,
    QSS_ref,
    LowQscore,
    // indels only:
    Repeat,
    iHpol,
    IndelBCNoise,
    QSI_ref,
    SIZE
};

inline
const char*
get_label(const unsigned idx)
{
    switch (idx)
    {
    case HighDepth:
        return "HighDepth";
        //return "DP"; // old strelka workflow name
    case BCNoise:
        return "BCNoise";
    case SpanDel:
        return "SpanDel";
    case QSS_ref:
        return "QSS_ref";
    case Repeat:
        return "Repeat";
    case iHpol:
        return "iHpol";
    case IndelBCNoise:
        return "BCNoise";
    case QSI_ref:
        return "QSI_ref";
    case LowQscore:
            return "LowQscore";
    default:
        assert(false && "Unknown vcf filter id");
        return nullptr;
    }
}
}



struct strelka_shared_modifiers
{
    strelka_shared_modifiers()
    {
        clear();
    }

    void
    set_filter(const STRELKA_VCF_FILTERS::index_t i)
    {
        filters.set(i);
    }

    void
    write_filters(std::ostream& os) const;

    void
    clear()
    {
        filters.reset();
    }

    std::bitset<STRELKA_VCF_FILTERS::SIZE> filters;
};


std::ostream& operator<<(std::ostream& os,const strelka_shared_modifiers& shmod);

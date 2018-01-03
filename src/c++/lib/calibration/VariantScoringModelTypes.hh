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
 *      Author: Morten Kallberg
 */

#pragma once

#include <cassert>


namespace SCORING_CALL_TYPE
{
enum index_t
{
    GERMLINE,
    RNA,
    SOMATIC,
    SIZE
};

inline
const char*
get_label(const index_t i)
{
    switch (i)
    {
    case GERMLINE:
        return "Germline";
    case RNA:
        return "RNAseq";
    case SOMATIC:
        return "Somatic";
    default:
        assert(false && "Unknown scoring model call type.");
        return nullptr;
    }
}
}


namespace SCORING_VARIANT_TYPE
{
enum index_t
{
    SNV,
    INDEL,
    SIZE
};

inline
const char*
get_label(const index_t i)
{
    switch (i)
    {
    case SNV:
        return "SNV";
    case INDEL:
        return "INDEL";
    default:
        assert(false && "Unknown scoring model variant type.");
        return nullptr;
    }
}
}

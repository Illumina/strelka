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

///
/// \author Chris Saunders
///

#include "htsapi/bam_header_util.hh"

#include "blt_util/thirdparty_push.h"

#include "boost/tokenizer.hpp"

#include "blt_util/thirdparty_pop.h"


bool
check_header_compatibility(
    const bam_header_t* h1,
    const bam_header_t* h2)
{
    if (h1->n_targets != h2->n_targets)
    {
        return false;
    }

    for (int32_t i(0); i<h1->n_targets; ++i)
    {
        if (h1->target_len[i] != h2->target_len[i]) return false;
        if (0 != strcmp(h1->target_name[i],h2->target_name[i])) return false;
    }
    return true;
}



std::string
get_bam_header_sample_name(
    const std::string& bam_header_text,
    const char* default_sample_name)
{
    using namespace boost;
    char_separator<char> sep("\t\n");
    tokenizer< char_separator<char>> tokens(bam_header_text, sep);
    for (const auto& t : tokens)
    {
        const auto sm = t.find("SM:");
        if (std::string::npos == sm) continue;
        return t.substr(sm+3);
    }

    return default_sample_name;
}

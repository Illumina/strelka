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

#include "htsapi/bam_header_util.hh"

#include "boost/tokenizer.hpp"



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

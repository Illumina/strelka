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

#include "ploidy_util.hh"



boost::optional<int>
parsePloidyFromBed(const char* line)
{
    boost::optional<int> result;

    if (line == nullptr) return result;

    unsigned tabcount(0);
    while (true)
    {
        if (*line=='\0' || *line=='\n') return result;
        if (*line=='\t') tabcount++;
        line++;
        if (tabcount>=4) break;
    }

    const char* s(line);
    int val = illumina::blt_util::parse_int(s);
    if (s != line) result.reset(val);

    return result;
}

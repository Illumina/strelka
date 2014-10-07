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

/// \file
///
/// \author Chris Saunders
///

#pragma once

#include <vector>


/// helps to control realignment complexity.
///
/// complexity is restricted as a function of the number of
/// indel candidates intersecting a read
///
struct starling_align_limit
{
    starling_align_limit(const unsigned max_alignment_count);

    unsigned
    get_max_toggle(const unsigned n_usable_indels) const
    {
        if (n_usable_indels>=_max_toggle.size())
        {
            return 1;
        }
        else
        {
            return _max_toggle[n_usable_indels];
        }
    }

private:
    std::vector<unsigned> _max_toggle;
};


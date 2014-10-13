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

#include <cstdint>


struct SiteNoise
{
    void
    clear()
    {
        total = 0;
        noise = 0;
        noise2 = 0;
    }

    double
    nfrac() const
    {
        return static_cast<double>(noise)/total;
    }

    double
    n2frac() const
    {
        return static_cast<double>(noise2)/total;
    }

    uint16_t total = 0;
    uint16_t noise = 0;
    uint16_t noise2 = 0;
};


void
set_noise_from_vcf(
    const char* line,
    SiteNoise& sn);


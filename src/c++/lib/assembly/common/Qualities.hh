// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Rumovsky
// Copyright (c) 2013-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//
/// \file
/// \author Ole Schulz-Trieglaff
///

#pragma once

static const int PHRED_OFFSET = 33;

// Phred score transformations
inline int char2phred(const char b)
{
    uint8_t v = b;
    assert(v >= PHRED_OFFSET);
    return v - PHRED_OFFSET;
}

inline char phred2char(const int p)
{
    uint8_t v = (p <= PHRED_OFFSET) ? p : 93;
    char c = v + PHRED_OFFSET;
    return c;
}


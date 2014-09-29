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

/// \brief provide access to cmake project version numbers

#pragma once

#include "common/config.h"

namespace starka
{

inline
const char*
getVersion()
{
    return STARKA_VERSION;
}


inline
const char*
getFullVersion()
{
    return STARKA_FULL_VERSION;
}

}

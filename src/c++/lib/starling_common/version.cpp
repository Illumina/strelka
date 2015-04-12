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

#include "version.hh"

#include "common/config.h"


const char*
getVersion()
{
    return WORKFLOW_VERSION;
}

inline
const char*
buildTime()
{
    return BUILD_TIME;
}

inline
const char*
cxxCompilerName()
{
    return CXX_COMPILER_NAME;
}

inline
const char*
compilerVersion()
{
    return COMPILER_VERSION;
}

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

#include "starling_common/prog_info_base.hh"
#include "starling_common/version.hh"

#include <string>


const char*
prog_info_base::
version() const
{
    static const std::string versionStr(starka::getVersion() + std::string(" ") + starka::buildTime());
    return versionStr.c_str();
}

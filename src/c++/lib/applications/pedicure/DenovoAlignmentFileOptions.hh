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

#include "PedicureSampleInfo.hh"

#include <string>
#include <vector>


struct  DenovoAlignmentFileOptions
{
    std::vector<std::string> alignmentFilename;

    // for now assume a 1-1 mapping from alignment files to samples:
    SampleInfoManager alignmentSampleInfo;
};


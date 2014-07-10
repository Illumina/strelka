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


namespace STAGE
{
enum index_t
{
    HEAD,
    READ_BUFFER,
    POST_ALIGN,
    POST_REGION, // haplotype specific stage
    POST_READ,  // haplotype specific stage
    POST_CALL,
    SIZE
};
}

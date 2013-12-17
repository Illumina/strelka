// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

/// \file

/// \author Chris Saunders
///
#ifndef __MATCH_DESCRIPTOR_HH
#define __MATCH_DESCRIPTOR_HH

/// \brief expand match descriptor to get reference sequence
///
/// ref is the corresponding reference base for each site in 'read'. A
/// gap character '-' is given for each insertion site in the read.
///
/// pos is the 1-indexed mapping offset from read to reference, in
/// read orientation. 0 indicates positions where the read is inserted
/// relative to reference
///
/// ref_read_length is the length on the reference sequence from the
/// first to the last mapped position of the read
///
void
expand_match_descriptor(const char* const read,
                        const char* const md,
                        char* ref,
                        unsigned* pos,
                        const unsigned read_length,
                        unsigned& ref_read_length);

#endif

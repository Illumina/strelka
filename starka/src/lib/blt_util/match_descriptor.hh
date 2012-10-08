// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2012 Illumina, Inc.
//
// This software is covered by the "Illumina Genome Analyzer Software
// License Agreement" and the "Illumina Source Code License Agreement",
// and certain third party copyright/licenses, and any user of this
// source file is bound by the terms therein (see accompanying files
// Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
// Illumina_Source_Code_License_Agreement.pdf and third party
// copyright/license notices).
//
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

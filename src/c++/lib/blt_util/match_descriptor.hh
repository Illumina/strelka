// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2016 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

/// \file

/// \author Chris Saunders
///

#pragma once

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

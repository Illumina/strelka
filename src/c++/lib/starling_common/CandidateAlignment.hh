//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2018 Illumina, Inc.
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

#pragma once

#include "alignment.hh"
#include "starling_common/indel_set.hh"
#include "starling_read_segment.hh"

#include <cassert>


/// \brief A rich alignment format used by the realignment routine
///
/// compared to a regular alignment this adds the insert sequence
/// for each indel (by keeping indel keys instead of a CIGAR), and
/// special fields for (possibly incomplete) edge indels.
///
struct CandidateAlignment
{
    bool
    operator<(const CandidateAlignment& rhs) const
    {
        if (al < rhs.al) return true;
        if (not (al == rhs.al)) return false;
        if (_indels < rhs._indels) return true;
        if (_indels != rhs._indels) return false;
        if (leading_indel_key < rhs.leading_indel_key) return true;
        if (not (leading_indel_key == rhs.leading_indel_key)) return false;
        return (trailing_indel_key < rhs.trailing_indel_key);
    }

    /// a new extension to specify the exact keys of all indels included
    /// in the alignment. this is a transitional bandaid required to
    /// distinguish insertions with different sequences, which can't
    /// be expressed in the path of the alignment structure below
    ///
    void
    setIndels(
        const indel_set_t& indels)
    {
        assert(not _isIndelsSet);
        _indels=indels;
        _isIndelsSet = true;
    }

    const indel_set_t&
    getIndels() const
    {
        assert(_isIndelsSet);
        return _indels;
    }

    alignment al;
    IndelKey leading_indel_key;
    IndelKey trailing_indel_key;

private:
    indel_set_t _indels;
    bool _isIndelsSet = false;
};


std::ostream& operator<<(std::ostream& os, const CandidateAlignment& cal);



/// \brief Get the keys of the indels present in the candidate alignment
///
void
getAlignmentIndels(
    const CandidateAlignment& cal,
    const reference_contig_segment& ref,
    const read_segment& rseg,
    const unsigned max_indel_size,
    const bool includeMismatches,
    indel_set_t& indels);

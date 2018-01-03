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

///
/// \author Chris Saunders
///

#pragma once


#include "blt_util/align_path.hh"
#include "blt_util/blt_types.hh"

#include <iosfwd>

//
// base class for non-chimeric alignments -- should apply to contigs
// and reads
//
// note read and alignment are always stored in fwd orientation:
//
struct alignment
{
    bool
    empty() const
    {
        return path.empty();
    }

    /// is there an indel exceeding max_indel_size?
    bool
    is_overmax(const unsigned max_indel_size) const;

    /// is there an adjacent insertion/deletion event in this
    /// alignment?
    bool
    is_seq_swap() const
    {
        return ALIGNPATH::is_seq_swap(path);
    }

    bool
    is_realignable(const unsigned max_indel_size) const
    {
        return (! is_overmax(max_indel_size));
    }

    void
    clear()
    {
        pos=0;
        path.clear();
        is_fwd_strand=true;
    }

    bool
    operator<(const alignment& rhs) const
    {
        if (pos<rhs.pos) return true;
        if (pos!=rhs.pos) return false;
        if (is_fwd_strand<rhs.is_fwd_strand) return true;
        if (is_fwd_strand!=rhs.is_fwd_strand) return false;

        const unsigned ps(path.size());
        const unsigned rps(rhs.path.size());
        if (ps<rps) return true;
        if (ps!=rps) return false;
        for (unsigned i(0); i<ps; ++i)
        {
            if (path[i]<rhs.path[i]) return true;
            if (not (path[i]==rhs.path[i])) return false;
        }
        return false;
    }

    bool
    operator==(const alignment& rhs) const
    {
        return ((pos==rhs.pos) && (path==rhs.path) && (is_fwd_strand==rhs.is_fwd_strand));
    }

    /////
    ALIGNPATH::path_t path;
    pos_t pos = 0;
    bool is_fwd_strand = true;
};


std::ostream& operator<<(std::ostream& os,const alignment& al);

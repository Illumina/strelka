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

    // is there an indel exceeding max_indel_size?
    bool
    is_overmax(const unsigned max_indel_size) const;

    // is there an adjacent insertion/deletion event in this
    // alignment?
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
        if (pos==rhs.pos)
        {
            if (is_fwd_strand<rhs.is_fwd_strand) return true;
            if (is_fwd_strand==rhs.is_fwd_strand)
            {
                const unsigned ps(path.size());
                const unsigned rps(rhs.path.size());
                if (ps<rps) return true;
                if (ps==rps)
                {
                    for (unsigned i(0); i<ps; ++i)
                    {
                        if (path[i]<rhs.path[i]) return true;
                        if (path[i]==rhs.path[i]) continue;
                        return false;
                    }
                }
            }
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

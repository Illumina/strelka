// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/downloads/sequencing/licenses/>.
//

/// \file
///
/// \author Chris Saunders
///

#ifndef __DEPTH_BUFFER_HH
#define __DEPTH_BUFFER_HH

#include "blt_util/blt_types.hh"

#include <cassert>

#include <map>


struct depth_buffer {

    unsigned
    val(const pos_t pos) const {
        const citer i(_data.find(pos));
        if(i == _data.end()) return 0;
        else                 return i->second;
    }

    void
    inc(const pos_t pos) {
        const iter i(_data.find(pos));
        if(i == _data.end()) _data[pos] = 1;
        else                 i->second += 1;
    }

    void
    clear_pos(const pos_t pos) {
        _data.erase(pos);
    }

    bool
    is_range_ge_than(const pos_t begin,
                     const pos_t end,
                     const unsigned depth) const {

        assert(begin <= end);
        citer i(_data.lower_bound(begin));
        const citer i_end(_data.upper_bound(end));
        for(; i!=i_end; ++i) if(i->second >= depth) return true;
        return false;
    }

private:
    typedef std::map<pos_t,unsigned> count_t;
    typedef count_t::iterator iter;
    typedef count_t::const_iterator citer;

    count_t _data;
};


#endif

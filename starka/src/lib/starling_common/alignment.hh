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
///
/// \author Chris Saunders
///

#ifndef __ALIGNMENT_HH
#define __ALIGNMENT_HH


#include "blt_util/blt_types.hh"
#include "starling_common/align_path.hh"

#include <iosfwd>

//
// base class for non-chimeric alignments -- should apply to contigs
// and reads
//
// note read and alignment are always stored in fwd orientation:
//
struct alignment {

    alignment() 
        : pos(0),
          is_fwd_strand(true) {}

    bool
    empty() const { return path.empty(); }

    // is there an (internal) indel exceeding MAX_INDEL_SIZE or
    // another reason this alignment can't be handled?
    bool
    is_overmax(const unsigned max_indel_size) const;

    // is there an adjacent insertion/deletion event in this
    // alignment?
    bool
    is_seq_swap() const { return ALIGNPATH::is_seq_swap(path); }

    bool
    is_realignable(const unsigned max_indel_size) const {
        return (! is_overmax(max_indel_size));
    }

    void
    clear() {
        pos=0;
        path.clear();
        is_fwd_strand=true;
    }

    bool
    operator<(const alignment& rhs) const {
        if(pos<rhs.pos) return true;
        if(pos==rhs.pos) {
            if(is_fwd_strand<rhs.is_fwd_strand) return true;
            if(is_fwd_strand==rhs.is_fwd_strand) {
                const unsigned ps(path.size());
                const unsigned rps(rhs.path.size());
                if(ps<rps) return true;
                if(ps==rps){
                    for(unsigned i(0);i<ps;++i){
                        if(path[i]<rhs.path[i]) return true;
                        if(path[i]==rhs.path[i]) continue;
                        return false;
                    }
                }
            }
        }
        return false;
    }

    bool
    operator==(const alignment& rhs) const {
        return ((pos==rhs.pos) and (path==rhs.path) and (is_fwd_strand==rhs.is_fwd_strand));
    }

    /////
    ALIGNPATH::path_t path;
    pos_t pos;
    bool is_fwd_strand;
};


std::ostream& operator<<(std::ostream& os,const alignment& al);


#endif

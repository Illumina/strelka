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
#ifndef __REFERENCE_CONTIG_SEGMENT_HH
#define __REFERENCE_CONTIG_SEGMENT_HH

#include "blt_util/blt_types.hh"

#include <string>


/// Manages a partial reference sequence segment
///
/// This object holds the reference sequence specified by the current
/// runs begin and end range, plus some padding on each side. To get 
/// this integrated into the current code as quickly as possible it
/// currently exposes the internal string object holding the sequence
/// data. When time allows this will be restricted so that a compressed
/// internal object can be used.
///
struct reference_contig_segment {

    reference_contig_segment()
        : _offset(0)
    {}

    char
    get_base(const pos_t pos) const {
        if(pos<_offset || pos>=end()) return 'N';
        return _seq[pos-_offset];
    }

    void
    get_substring(const pos_t pos,
                  const pos_t length,
                  std::string& substr) const {

        if(pos<_offset || (pos+length)>end()) {
            //slow path (minority of calls):
            substr.clear();
            if(pos<_offset) {
                substr += std::string(std::min(length,(_offset-pos)),'N');
            }
            
            const unsigned rstart(std::min(pos,_offset));
            const unsigned rend(std::max(pos+length,end()));
            if(rend>rstart) {
                substr += _seq.substr(rstart,rend-rstart);
            }
            if((pos+length)>end()) {
                substr += std::string(std::min(length,((pos+length)-end())),'N');
            }
        } else {
            //fast path
            substr.assign(_seq,pos-_offset,length);
        }
    }

    std::string& seq() { return _seq; }
    const std::string& seq() const { return _seq; }

    void
    set_offset(const pos_t offset) {
        _offset=offset;
    }

    pos_t
    end() const { return _offset+_seq.size(); }

private:

    pos_t _offset;
    std::string _seq;
};


#endif

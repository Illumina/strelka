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


#include "blt_util/read_util.hh"

#include <cstring>



void
get_read_align_strand_end_skip(const char* const read,
                               const unsigned read_size,
                               unsigned& end_skip) {

    unsigned read_end(read_size);

    while (read_end>0) {
        if (read[read_end-1]=='N') read_end--;
        else break;
    }

    end_skip=read_size-read_end;
}



void
get_read_fwd_strand_skip(const char* const read,
                         const unsigned read_size,
                         const bool is_fwd_strand,
                         unsigned& begin_skip,
                         unsigned& end_skip) {

    begin_skip=0;
    if (is_fwd_strand) {
        get_read_align_strand_end_skip(read,read_size,end_skip);
    } else {
        end_skip=0;
        while (begin_skip<read_size) {
            if (read[begin_skip]=='N') begin_skip++;
            else break;
        }
    }
}



static
void
get_read_align_strand_end_skip(const bam_seq& bseq,
                               unsigned& end_skip) {

    unsigned read_end(bseq.size());

    while (read_end>0) {
        if (bseq.get_char(read_end-1)=='N') read_end--;
        else break;
    }

    end_skip=bseq.size()-read_end;
}



void
get_read_fwd_strand_skip(const bam_seq& bseq,
                         const bool is_fwd_strand,
                         unsigned& begin_skip,
                         unsigned& end_skip) {

    begin_skip=0;
    if (is_fwd_strand) {
        get_read_align_strand_end_skip(bseq,end_skip);
    } else {
        end_skip=0;
        const unsigned bsize(bseq.size());
        while (begin_skip<bsize) {
            if (bseq.get_char(begin_skip)=='N') begin_skip++;
            else break;
        }
    }
}

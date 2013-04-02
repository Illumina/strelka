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

/// \author Chris Saunders
///

#include "blt_util/bam_streamer.hh"
#include "blt_util/log.hh"

#include <cstdlib>

#include <iostream>



bam_streamer::
bam_streamer(const char* filename,
             const char* region)
    : _is_record_set(false), _record_no(0), _stream_name(filename), _is_region(false),
      _bfp(NULL), _bidx(NULL), _biter(NULL) {

    assert(NULL != filename);
    assert('\0' != *filename);

    _bfp = samopen(filename, "rb", 0);

    if(NULL == _bfp) {
        log_os << "ERROR: Failed to open SAM/BAM file: " << filename << "\n";
        exit(EXIT_FAILURE);
    }

    // read the whole BAM file:
    if(NULL == region) return;

    _is_region=true;
    _region=region;

    // use the BAM index to read a region of the BAM file
    if(! (_bfp->type&0x01)) {
        log_os << "ERROR: file must be in BAM format for region lookup: " << filename << "\n";
        exit(EXIT_FAILURE);
    }

    /// TODO: Find out whether _bidx can be destroyed after the BAM
    /// iterator is created, in which case this could be a local
    /// variable. Until we know, _bidx should persist for the lifetime
    /// of _biter
    _bidx = bam_index_load(filename); // load BAM index
    if (NULL == _bidx) {
        log_os << "ERROR: BAM index is not available for file: " << filename << "\n";
        exit(EXIT_FAILURE);
    }

    int ref,beg,end;
    bam_parse_region(_bfp->header, region, &ref, &beg, &end); // parse the region

    if (ref < 0) {
        log_os << "ERROR: Invalid region: '" <<  region << "' specified for BAM file: " << filename << "\n";
        exit(EXIT_FAILURE);
    }

    _biter = bam_iter_query(_bidx,ref,beg,end);
}



bam_streamer::
~bam_streamer() {
    if(NULL != _biter) bam_iter_destroy(_biter);
    if(NULL != _bidx) bam_index_destroy(_bidx);
    if(NULL != _bfp) samclose(_bfp);
}



bool
bam_streamer::
next() {
    if(NULL==_bfp) return false;

    int ret;
    if(NULL == _biter) {
        ret = samread(_bfp, _brec._bp);
    } else {
        ret = bam_iter_read(_bfp->x.bam, _biter, _brec._bp);
    }

    _is_record_set=(ret >= 0);
    if(_is_record_set) _record_no++;

    return _is_record_set;
}



const char*
bam_streamer::
target_id_to_name(const int32_t tid) const {
    // assert(tid < _bfp->header->n_targets);
    if(tid<0) {
        static const char unmapped[] = "*";
        return unmapped;
    }
    return _bfp->header->target_name[tid];
}



int32_t
bam_streamer::
target_name_to_id(const char* seq_name) const {
    return bam_get_tid(_bfp->header,seq_name);
}



void
bam_streamer::
report_state(std::ostream& os) const {

    const bam_record* bamp(get_record_ptr());

    os << "\tbam_stream_label: " << name() << "\n";
    if(_is_region) {
        os << "\tbam_stream_selected_region: " << _region << "\n";
    }
    if(NULL != bamp) {
        os << "\tbam_stream_record_no: " << record_no() << "\n";
        os << "\tbam_record QNAME/read_number: " << bamp->qname() << "/" << bamp->read_no() << "\n";
        const char* chrom_name(target_id_to_name(bamp->target_id()));
        os << "\tbam record RNAME: " << chrom_name << "\n";
        os << "\tbam record POS: " << bamp->pos() << "\n";
    } else {
        os << "\tno bam record currently set\n";
    }
}

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
/// object which accepts as input bam and contig files from multiple
/// samples and presents them in the order expected by
/// starling_pos_processor
///
///
/// \author Chris Saunders
///
/// note coding convention for all ranges '_pos fields' is:
/// XXX_begin_pos is zero-indexed position at the begining of the range
/// XXX_end_pos is zero-index position 1 step after the end of the range
///

#pragma once


#include "blt_util/id_map.hh"
#include "blt_util/bam_streamer.hh"
#include "blt_util/vcf_streamer.hh"
#include "starling_common/starling_types.hh"
#include "starling_common/grouper_contig_util.hh"

#include <map>
#include <queue>


namespace INPUT_TYPE {
enum index_t { NONE, READ, CONTIG, INDEL, FORCED_OUTPUT };
}

struct starling_input_stream_hander;

struct starling_input_stream_data {

    void
    register_reads(bam_streamer& bs,
                   const sample_id_t sample_no = 0) {
        if (_reads.test_key(sample_no)) register_error("reads",sample_no);
        _reads.insert(sample_no,&bs);
    }

    void
    register_contigs(contig_reader& cr,
                     const sample_id_t sample_no = 0) {
        if (_contigs.test_key(sample_no)) register_error("contigs",sample_no);
        _contigs.insert(sample_no,&cr);
    }

    void
    register_indels(vcf_streamer& vr,
                    const sample_id_t sample_no = 0) {
        // unlike reads/contigs, we allow multiple files associated with the same
        // sample_no for input indels:
        _indels.push_back(std::make_pair(sample_no,&vr));
    }

    void
    register_forced_output(vcf_streamer& vr,
                           const sample_id_t sample_no = 0) {
        // sites and indels in these files must be included in the snv/indel output, this means that
        // any indels in these files are also candidate indels:
        _output.push_back(std::make_pair(sample_no,&vr));
    }

private:

    void
    register_error(const char* label,
                   const sample_id_t sample_no) const;


/////////// data:
    friend struct starling_input_stream_handler;
    typedef id_map<sample_id_t,bam_streamer*> reads_t;
    typedef id_map<sample_id_t,contig_reader*> contigs_t;
    typedef std::vector<std::pair<sample_id_t, vcf_streamer*> > indels_t;

    reads_t _reads;
    contigs_t _contigs;
    indels_t _indels;
    indels_t _output;
};



/// abstracts different record types (bam/contig/vcf, etc...) so that these can be
/// sorted and handled in order
///
struct input_record_info {

    input_record_info(const pos_t p = 0,
                      const INPUT_TYPE::index_t t = INPUT_TYPE::NONE,
                      const sample_id_t i = 0,
                      const unsigned s = 0)
        :  pos(p), itype(t), sample_no(i), _order(s) {}

    // reverse logic implied by operator< such that the 'lower' values
    // we'd like to see first will come up on top of the
    // priority_queue
    //
    bool
    operator<(const input_record_info& rhs) const {
        if (pos > rhs.pos) return true;
        if (pos == rhs.pos) {
            if (itype < rhs.itype) return true;
            if (itype==rhs.itype) {
                if (sample_no > rhs.sample_no) return true;
                if (sample_no == rhs.sample_no) {
                    return (_order > rhs._order);
                }
            }
        }
        return false;
    }

    unsigned get_order() const { return _order; }

    pos_t pos;
    INPUT_TYPE::index_t itype;
    sample_id_t sample_no;

private:
    friend struct starling_input_stream_handler;

    // record the submission order:
    unsigned _order;
};



// streams multiple bams, contig and vcf files to present the data
// in positional order (but with offsets for contigs and vcfs to
// run ahead of the bam reads)
//
struct starling_input_stream_handler {

    starling_input_stream_handler(const starling_input_stream_data& data,
                                  const pos_t contig_lead = 1000,
                                  const pos_t indel_lead = 100,
                                  const pos_t output_lead = 100);

    bool next();

    input_record_info
    get_current() const { return _current; }

    pos_t
    get_head_pos() const { return _head_pos; }

private:

    void
    push_next(const INPUT_TYPE::index_t itype,
              const sample_id_t sample_no,
              const unsigned order);


///////////////////////////////// data:
    const starling_input_stream_data& _data;

    // contig_lead controls the amount by which we read the
    // local-assembly contig and contig read buffer ahead of the
    // genomic reads:
    //
    // \TODO redesign so that we handle the much larger deletions coming
    // from GROUPER (more) correctly. Idea is (1) to prevent any possible
    // double counting of a widely separated GROUPER and genomic read and
    // (2) to handle deletion spanning reads as segmented reads using the
    // exon handling system.
    //
    const pos_t _contig_lead;

    // ditto for indel_lead with indels and forced output:
    const pos_t _indel_lead;
    const pos_t _output_lead;

    input_record_info _current;
    input_record_info _last;

    bool _is_end;

    bool _is_head_pos;
    pos_t _head_pos;

    std::priority_queue<input_record_info> _stream_queue;
};


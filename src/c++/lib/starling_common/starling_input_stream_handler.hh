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
/// note coding convention for all ranges '_pos fields' is:
/// XXX_begin_pos is zero-indexed position at the beginning of the range
/// XXX_end_pos is zero-index position 1 step after the end of the range
///

#pragma once


#include "blt_util/id_map.hh"
#include "blt_util/bam_streamer.hh"
#include "blt_util/vcf_streamer.hh"
#include "starling_common/starling_types.hh"

#include <map>
#include <queue>


namespace INPUT_TYPE
{
enum index_t {
    NONE,
    READ,
    INDEL,
    FORCED_OUTPUT,
    NOISE
};
}

struct starling_input_stream_hander;

struct starling_input_stream_data
{
    void
    register_reads(bam_streamer& bs,
                   const sample_id_t sample_no = 0)
    {
        if (_reads.test_key(sample_no)) register_error("reads",sample_no);
        _reads.insert(sample_no,&bs);
    }

    /// unlike reads/contigs, we allow multiple files associated with the same
    /// sample_no for input indels:
    void
    register_indels(vcf_streamer& vr,
                    const sample_id_t sample_no = 0)
    {
        _indels.push_back(std::make_pair(sample_no,&vr));
    }

    /// sites and indels in these files must be included in the snv/indel output, this means that
    /// any indels in these files are also candidate indels:
    void
    register_forced_output(vcf_streamer& vr,
                           const sample_id_t sample_no = 0)
    {
        _output.push_back(std::make_pair(sample_no,&vr));
    }

    /// sites and indels in these files will be used to estimate low-freqeuncy noise
    void
    register_noise(
        vcf_streamer& vr,
        const sample_id_t sample_no = 0)
    {
        _noise.push_back(std::make_pair(sample_no,&vr));
    }

private:

    void
    register_error(const char* label,
                   const sample_id_t sample_no) const;


/////////// data:
    friend struct starling_input_stream_handler;
    typedef id_map<sample_id_t,bam_streamer*> reads_t;
    typedef std::vector<std::pair<sample_id_t, vcf_streamer*> > indels_t;

    reads_t _reads;
    indels_t _indels;
    indels_t _output;
    indels_t _noise;
};



/// abstracts different record types (bam/contig/vcf, etc...) so that these can be
/// sorted and handled in order
///
struct input_record_info
{
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
    operator<(const input_record_info& rhs) const
    {
        if (pos > rhs.pos) return true;
        if (pos == rhs.pos)
        {
            if (itype < rhs.itype) return true;
            if (itype==rhs.itype)
            {
                if (sample_no > rhs.sample_no) return true;
                if (sample_no == rhs.sample_no)
                {
                    return (_order > rhs._order);
                }
            }
        }
        return false;
    }

    unsigned get_order() const
    {
        return _order;
    }

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
struct starling_input_stream_handler
{
    starling_input_stream_handler(const starling_input_stream_data& data);

    bool next();

    input_record_info
    get_current() const
    {
        return _current;
    }

    pos_t
    get_head_pos() const
    {
        return _head_pos;
    }

private:

    void
    push_next(const INPUT_TYPE::index_t itype,
              const sample_id_t sample_no,
              const unsigned order);

///////////////////////////////// data:
    const starling_input_stream_data& _data;

    // vcf_lead controls the amount by which we read the
    // vcf buffer(s) ahead of the bam reads:
    //
    const pos_t _vcf_lead = 100;

    input_record_info _current;
    input_record_info _last;

    bool _is_end = false;

    bool _is_head_pos = false;
    pos_t _head_pos = 0;

    std::priority_queue<input_record_info> _stream_queue;
};


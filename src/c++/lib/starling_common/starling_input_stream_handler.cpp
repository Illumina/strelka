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

/// \file
///
/// \author Chris Saunders
///
/// note coding convention for all ranges '_pos fields' is:
/// XXX_begin_pos is zero-indexed position at the begining of the range
/// XXX_end_pos is zero-index position 1 step after the end of the range
///


#include "blt_util/blt_exception.hh"
#include "blt_util/log.hh"
#include "starling_common/starling_input_stream_handler.hh"

#include <cstdlib>

#include <iostream>
#include <sstream>



static
const char*
input_type_label(const INPUT_TYPE::index_t i)
{
    using namespace INPUT_TYPE;

    switch (i)
    {
    case NONE   :
        return "NONE";
    case READ   :
        return "READ";
    case INDEL  :
        return "INDEL";
    case FORCED_OUTPUT  :
        return "FORCED_OUTPUT";
    case NOISE  :
        return "NOISE";
    default :
        log_os << "ERROR: unrecognized event type.\n";
        exit(EXIT_FAILURE);
    }
}



void
starling_input_stream_data::
register_error(const char* label,
               const sample_id_t sample_no) const
{
    log_os << "ERROR: attempting to register " << label
           << " with sample number: " << sample_no
           << " more than once\n";
    exit(EXIT_FAILURE);
}




starling_input_stream_handler::
starling_input_stream_handler(const starling_input_stream_data& data)
    : _data(data)
{
    // initial loading for _stream_queue:
    const unsigned rs(_data._reads.size());
    for (unsigned i(0); i<rs; ++i)
    {
        push_next(INPUT_TYPE::READ,_data._reads.get_key(i),i);
    }
    const unsigned is(_data._indels.size());
    for (unsigned i(0); i<is; ++i)
    {
        push_next(INPUT_TYPE::INDEL,_data._indels[i].first,i);
    }
    const unsigned os(_data._output.size());
    for (unsigned i(0); i<os; ++i)
    {
        push_next(INPUT_TYPE::FORCED_OUTPUT,_data._output[i].first,i);
    }
    const unsigned ns(_data._noise.size());
    for (unsigned i(0); i<ns; ++i)
    {
        push_next(INPUT_TYPE::NOISE,_data._noise[i].first,i);
    }
}



bool
starling_input_stream_handler::
next()
{
    if (_is_end) return false;

    while (true)
    {
        if (_current.itype != INPUT_TYPE::NONE)
        {
            // reload stream_queue with current type and sample_no;
            push_next(_current.itype,_current.sample_no,_current._order);
            _last=_current;
        }

        if (_stream_queue.empty())
        {
            _current=input_record_info();
            _is_end=true;
            return false;
        }
        bool is_usable(true);
        _current=_stream_queue.top();
        _stream_queue.pop();

        if (_is_head_pos &&
            (_current.pos < _head_pos))
        {
            if (_current.itype == INPUT_TYPE::READ)
            {
                std::ostringstream oss;
                oss << "ERROR: unexpected read order:\n"
                    << "\tInput-record with pos/type/sample_no: "
                    << (_current.pos+1) << "/" << input_type_label(_current.itype) << "/" << _current.sample_no
                    << " follows pos/type/sample_no: "
                    << (_last.pos+1) << "/" << input_type_label(_last.itype) << "/" << _current.sample_no << "\n";
                throw blt_exception(oss.str().c_str());
            }
            else if ((_current.itype == INPUT_TYPE::INDEL) ||
                     (_current.itype == INPUT_TYPE::FORCED_OUTPUT) ||
                     (_current.itype == INPUT_TYPE::NOISE))
            {
                std::ostringstream oss;
                oss << "ERROR: unexpected vcf record order:\n"
                    << "\tInput-record with pos/type/sample_no: "
                    << (_current.pos+1) << "/" << input_type_label(_current.itype) << "/" << _current.sample_no
                    << " follows pos/type/sample_no: "
                    << (_last.pos+1) << "/" << input_type_label(_last.itype) << "/" << _current.sample_no << "\n";
                throw blt_exception(oss.str().c_str());
            }
            else
            {
                assert(false && "Unknown input type");
            }
        }

        if (_is_head_pos)
        {
            _head_pos=std::max(_head_pos,_current.pos);
        }
        else
        {
            _is_head_pos=true;
            _head_pos=_current.pos;
        }

        if (is_usable) break;
    }
    return true;
}



static
void
get_next_read_pos(bool& is_next_read,
                  pos_t& next_read_pos,
                  bam_streamer& read_stream)
{
    is_next_read=read_stream.next();
    if (is_next_read)
    {
        const bam_record& read_rec(*(read_stream.get_record_ptr()));
        next_read_pos=(read_rec.pos()-1);
    }
    else
    {
        next_read_pos=0;
    }
}



//
static
void
get_next_indel_pos(bool& is_next_indel,
                   pos_t& next_indel_pos,
                   vcf_streamer& indel_stream)
{
    static const bool is_indel_only(true);
    is_next_indel=indel_stream.next(is_indel_only);
    if (is_next_indel)
    {
        const vcf_record& vcf_rec(*(indel_stream.get_record_ptr()));
        next_indel_pos=(vcf_rec.pos-1);
    }
    else
    {
        next_indel_pos=0;
    }
}



//
static
void
get_next_forced_output_pos(bool& is_next_variant,
                           pos_t& next_variant_pos,
                           vcf_streamer& variant_stream)
{
    static const bool is_indel_only(false);
    is_next_variant=variant_stream.next(is_indel_only);
    if (is_next_variant)
    {
        const vcf_record& vcf_rec(*(variant_stream.get_record_ptr()));
        next_variant_pos=(vcf_rec.pos-1);
    }
    else
    {
        next_variant_pos=0;
    }
}



//
static
void
get_next_noise_pos(
    bool& is_next_variant,
    pos_t& next_variant_pos,
    vcf_streamer& variant_stream)
{
    static const bool is_indel_only(false);
    is_next_variant=variant_stream.next(is_indel_only);
    if (is_next_variant)
    {
        const vcf_record& vcf_rec(*(variant_stream.get_record_ptr()));
        next_variant_pos=(vcf_rec.pos-1);
    }
    else
    {
        next_variant_pos=0;
    }
}



void
starling_input_stream_handler::
push_next(const INPUT_TYPE::index_t itype,
          const sample_id_t sample_no,
          const unsigned order)
{
    bool is_next(false);
    pos_t next_pos;
    if       (itype == INPUT_TYPE::READ)
    {
        bam_streamer& read_stream(*(_data._reads.get_value(order)));
        get_next_read_pos(is_next,next_pos,read_stream);
    }
    else
    {
        if (itype == INPUT_TYPE::INDEL)
        {
            vcf_streamer& vcf_stream(*(_data._indels[order].second));
            get_next_indel_pos(is_next,next_pos,vcf_stream);
        }
        else if (itype == INPUT_TYPE::FORCED_OUTPUT)
        {
            vcf_streamer& vcf_stream(*(_data._output[order].second));
            get_next_forced_output_pos(is_next,next_pos,vcf_stream);
        }
        else if (itype == INPUT_TYPE::NOISE)
        {
            vcf_streamer& vcf_stream(*(_data._noise[order].second));
            get_next_noise_pos(is_next,next_pos,vcf_stream);
        }
        else
        {
            assert(false && "Unknown input type");
        }
        next_pos -= std::min(_vcf_lead,next_pos);
    }
    if (! is_next) return;
    _stream_queue.push(input_record_info(next_pos,itype,sample_no,order));
}


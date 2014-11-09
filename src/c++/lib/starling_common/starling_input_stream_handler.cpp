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
input_type_label(
    const INPUT_TYPE::index_t i)
{
    using namespace INPUT_TYPE;

    switch (i)
    {
    case NONE :
        return "NONE";
    case READ :
        return "READ";
    case INDEL :
        return "INDEL";
    case FORCED_OUTPUT :
        return "FORCED_OUTPUT";
    case PLOIDY_REGION :
        return "PLOIDY_REGION";
    case NOCOMPRESS_REGION :
        return "NOCOMPRESS_REGION";
    case NOISE :
        return "NOISE";
    default :
        log_os << "ERROR: unrecognized event type.\n";
        exit(EXIT_FAILURE);
    }
}



void
starling_input_stream_data::
register_error(
    const char* label,
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
    const unsigned ps(_data._ploidy.size());
    for (unsigned i(0); i<ps; ++i)
    {
        push_next(INPUT_TYPE::PLOIDY_REGION,_data._ploidy[i].first,i);
    }
    const unsigned cs(_data._nocompress.size());
    for (unsigned i(0); i<cs; ++i)
    {
        push_next(INPUT_TYPE::NOCOMPRESS_REGION,_data._nocompress[i].first,i);
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
            std::ostringstream oss;
            if (_current.itype == INPUT_TYPE::READ)
            {
                oss << "ERROR: unexpected read order:\n"
                    << "\tInput-record with pos/type/sample_no: "
                    << (_current.pos+1) << "/" << input_type_label(_current.itype) << "/" << _current.sample_no
                    << " follows pos/type/sample_no: "
                    << (_last.pos+1) << "/" << input_type_label(_last.itype) << "/" << _current.sample_no << "\n";
            }
            else if ((_current.itype == INPUT_TYPE::INDEL) ||
                     (_current.itype == INPUT_TYPE::FORCED_OUTPUT) ||
                     (_current.itype == INPUT_TYPE::NOISE))
            {
                oss << "ERROR: unexpected vcf record order:\n"
                    << "\tInput-record with pos/type/sample_no: "
                    << (_current.pos+1) << "/" << input_type_label(_current.itype) << "/" << _current.sample_no
                    << " follows pos/type/sample_no: "
                    << (_last.pos+1) << "/" << input_type_label(_last.itype) << "/" << _current.sample_no << "\n";
            }
            else if ((_current.itype == INPUT_TYPE::PLOIDY_REGION) ||
                     (_current.itype == INPUT_TYPE::NOCOMPRESS_REGION))
            {
                oss << "ERROR: unexpected bed record order:\n"
                    << "\tInput-record with begin/type/sample_no: "
                    << (_current.pos+1) << "/" << input_type_label(_current.itype) << "/" << _current.sample_no
                    << " follows record with begin/type/sample_no: "
                    << (_last.pos+1) << "/" << input_type_label(_last.itype) << "/" << _current.sample_no << "\n";
            }
            else
            {
                assert(false && "Unknown input type");
            }
            throw blt_exception(oss.str().c_str());
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
get_next_variant_pos(
    const bool is_indel_only,
    bool& is_next_variant,
    pos_t& next_variant_pos,
    vcf_streamer& variant_stream)
{
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
get_next_region(
    bool& is_next_region,
    pos_t& next_region_pos,
    bed_streamer& region_stream)
{
    is_next_region=region_stream.next();
    if (is_next_region)
    {
        const bed_record& bed_rec(*(region_stream.get_record_ptr()));
        next_region_pos=(bed_rec.begin-1);
    }
    else
    {
        next_region_pos=0;
    }
}



void
starling_input_stream_handler::
push_next(
    const INPUT_TYPE::index_t itype,
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
    else if ((itype == INPUT_TYPE::PLOIDY_REGION) ||
             (itype == INPUT_TYPE::NOCOMPRESS_REGION))
    {
        bed_streamer* bed_stream(nullptr);
        if (itype == INPUT_TYPE::PLOIDY_REGION)
        {
            bed_stream=(_data._ploidy[order].second);
        }
        else if (itype == INPUT_TYPE::NOCOMPRESS_REGION)
        {
            bed_stream=(_data._nocompress[order].second);
        }
        else
        {
            assert(false && "Unknown input region type");
        }
        get_next_region(is_next,next_pos,*bed_stream);
        next_pos -= std::min(_bed_lead,next_pos);
    }
    else
    {
        vcf_streamer* vcf_stream(nullptr);
        bool is_indel_only(false);
        if (itype == INPUT_TYPE::INDEL)
        {
            vcf_stream=((_data._indels[order].second));
            is_indel_only=true;
        }
        else if (itype == INPUT_TYPE::FORCED_OUTPUT)
        {
            vcf_stream=((_data._output[order].second));
        }
        else if (itype == INPUT_TYPE::NOISE)
        {
            vcf_stream=((_data._noise[order].second));
        }
        else
        {
            assert(false && "Unknown input variant type");
        }
        get_next_variant_pos(is_indel_only,is_next,next_pos,*vcf_stream);
        next_pos -= std::min(_vcf_lead,next_pos);
    }
    if (! is_next) return;
    _stream_queue.push(input_record_info(next_pos,itype,sample_no,order));
}


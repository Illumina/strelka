/*
 * bedstreamprocessor.cpp
 *
 *  Created on: Jun 1, 2015
 *      Author: jduddy
 */

#include "bedstreamprocessor.hh"


bed_stream_processor::bed_stream_processor(const std::string& bed_file_name, const std::string& region)
: _bed_streamer(bed_file_name.c_str(), region.c_str())
, _region(region)
, _is_end(false)
{
    load_next_region();
}

void bed_stream_processor::load_next_region()
{
    if (!_is_end)
    {
        if (_bed_streamer.next())
        {
            _current_record = _bed_streamer.get_record_ptr();
            if (_current_record->chrom != _region)
            {
                _is_end = true;
                _current_record = nullptr;
            }
        }
        else
        {
            _current_record = nullptr;
            _is_end = true;
        }
    }
}

void bed_stream_processor::load_next_region_if_needed(pos_t position)
{
    while (_current_record != nullptr && _current_record->end < position)
        load_next_region();
}

bool bed_stream_processor::process(site_info& si)
{
    load_next_region_if_needed(si.pos);
    if (_current_record == nullptr || si.pos < _current_record->begin || si.pos >= _current_record->end)
    {
        // TODO: make sure we don't apply to stuff we shouldn't? No coverage? I dunno
        si.smod.set_filter(VCF_FILTERS::OffTarget);
    }
    return true;
}

bool bed_stream_processor::process(indel_info& ii)
{
    load_next_region_if_needed(ii.pos);
    // include indels if they start in a targeted region
    if (_current_record == nullptr || ii.pos < _current_record->begin || ii.pos >= _current_record->end)
    {
        ii.imod.set_filter(VCF_FILTERS::OffTarget);
    }
    return true;
}

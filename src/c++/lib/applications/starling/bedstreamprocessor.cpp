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
/*
 * bedstreamprocessor.cpp
 *
 *  Created on: Jun 1, 2015
 *      Author: jduddy
 */

#include "bedstreamprocessor.hh"


bed_stream_processor::bed_stream_processor(const std::string& bed_file_name, const std::string& region, std::shared_ptr<variant_pipe_stage_base> next_stage)
    : variant_pipe_stage_base(next_stage)
    , _region(region)
    , _is_end(false)
{
    if (!bed_file_name.empty())
    {
        _bed_streamer.reset(new bed_streamer(bed_file_name.c_str(), region.c_str()));
        load_next_region();
    }
}

void bed_stream_processor::load_next_region()
{
    if (!_is_end)
    {
        if (_bed_streamer->next())
        {
            _current_record = _bed_streamer->get_record_ptr();
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

void bed_stream_processor::process(std::unique_ptr<site_info> si)
{
    if (this->_bed_streamer)
    {
        load_next_region_if_needed(si->pos);
        if (_current_record == nullptr || si->pos < _current_record->begin || si->pos >= _current_record->end)
        {
            // TODO: make sure we don't apply to stuff we shouldn't? No coverage? I dunno
            si->set_filter(VCF_FILTERS::OffTarget);
        }
    }
    _sink->process(std::move(si));
}

void bed_stream_processor::process(std::unique_ptr<indel_info> ii)
{
    if (this->_bed_streamer)
    {
        load_next_region_if_needed(ii->pos);
        // include indels if they start in a targeted region
        if (_current_record == nullptr || ii->pos < _current_record->begin || ii->pos >= _current_record->end)
        {
            ii->set_filter(VCF_FILTERS::OffTarget);
        }
    }
    _sink->process(std::move(ii));
}

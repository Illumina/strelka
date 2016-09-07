// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2016 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//
/*
 * bedstreamprocessor.cpp
 *
 *  Created on: Jun 1, 2015
 *      Author: jduddy
 */

#include "bedstreamprocessor.hh"


bed_stream_processor::
bed_stream_processor(
    const std::string& bed_file_name,
    const std::string& region,
    std::shared_ptr<variant_pipe_stage_base> next_stage)
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



void
bed_stream_processor::
load_next_region()
{
    if (_is_end) return;

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



void
bed_stream_processor::
load_next_region_if_needed(const pos_t pos)
{
    while ((_current_record != nullptr) && (_current_record->end < pos))
    {
        load_next_region();
    }
}



void
bed_stream_processor::
processLocus(LocusInfo& locus)
{
    if (not _bed_streamer) return;

    load_next_region_if_needed(locus.pos);

    // include alleles if they start in a targeted region
    if (_current_record == nullptr || locus.pos < _current_record->begin || locus.pos >= _current_record->end)
    {
        locus.filters.set(GERMLINE_VARIANT_VCF_FILTERS::OffTarget);
    }
}



void
bed_stream_processor::
process(std::unique_ptr<GermlineSiteLocusInfo> si)
{
    processLocus(*si);
    _sink->process(std::move(si));
}



void
bed_stream_processor::
process(std::unique_ptr<GermlineIndelLocusInfo> ii)
{
    processLocus(*ii);
    _sink->process(std::move(ii));
}

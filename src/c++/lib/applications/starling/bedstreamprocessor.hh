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
 * bedstreamprocessor.hh
 *
 *  Created on: Jun 1, 2015
 *      Author: jduddy
 */

#pragma once
#include "htsapi/bed_streamer.hh"
#include "variant_pipe_stage_base.hh"

class bed_stream_processor : public variant_pipe_stage_base
{
public:
    bed_stream_processor(const std::string& bed_file_name, const std::string& region, std::shared_ptr<variant_pipe_stage_base> next_stage);

    void process(std::unique_ptr<site_info> si) override;
    void process(std::unique_ptr<indel_info> ii) override;
private:
    std::unique_ptr<bed_streamer> _bed_streamer;
    const bed_record* _current_record;

    // TODO: move enforcement of regions into bed_streamer
    std::string _region;
    bool _is_end;

    void load_next_region();
    void load_next_region_if_needed(pos_t position);
};


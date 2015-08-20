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
    bed_stream_processor(const std::string& bed_file_name, const std::string& region, variant_pipe_stage_base& next_stage);

    void process(site_info& si) override;
    void process(indel_info& ii) override;
private:
    std::unique_ptr<bed_streamer> _bed_streamer;
    const bed_record* _current_record;

    // TODO: move enforcement of regions into bed_streamer
    std::string _region;
    bool _is_end;

    void load_next_region();
    void load_next_region_if_needed(pos_t position);
};


/*
 * bedstreamprocessor.hh
 *
 *  Created on: Jun 1, 2015
 *      Author: jduddy
 */

#pragma once
#include "htsapi/bed_streamer.hh"
#include "variant_pipe_stage.hh"

class bed_stream_processor : public variant_pipe_stage
{
public:
    bed_stream_processor(const std::string& bed_file_name, const std::string& region, variant_pipe_stage& next_stage);

    void process(site_info& si) override;
    void process(indel_info& ii) override;
private:
    bed_streamer _bed_streamer;
    const bed_record* _current_record;

    // TODO: move enforcement of regions into bed_streamer
    std::string _region;
    bool _is_end;

    void load_next_region();
    void load_next_region_if_needed(pos_t position);
};


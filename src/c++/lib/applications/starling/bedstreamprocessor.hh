/*
 * bedstreamprocessor.hh
 *
 *  Created on: Jun 1, 2015
 *      Author: jduddy
 */

#pragma once
#include "variant_processor.hh"
#include "htsapi/bed_streamer.hh"

class bed_stream_processor : public variant_processor
{
public:
    bed_stream_processor(const std::string& bed_file_name, const std::string& region);

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


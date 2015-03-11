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
 * predictor.hh
 *
 *  Created on: Aug 10, 2013
 *  Author: Morten Kallberg
 */

#pragma once

#include "gvcf_locus_info.hh"
#include "starling_common/starling_base_shared.hh"
#include "starling_common/starling_read_buffer.hh"
#include "starling_common/starling_pos_processor_util.hh"
#include "blt_util/RegionTracker.hh"
#include "htsapi/bed_streamer.hh"
#include "applications/starling/starling_shared.hh"

#include <boost/algorithm/string.hpp>

//#define DEBUG_predictor


#ifdef DEBUG_predictor
#include "blt_util/log.hh"
#endif

struct predictor
{
    predictor(
              const starling_base_options& init_opt,
              const starling_deriv_options& init_dopt,
              const RegionTracker& init_nocompress_regions)
      : regions_file(init_opt.assembly_regions_filename), opt(init_opt), dopt(init_dopt),
        _nocompress_regions(init_nocompress_regions)
    {
        //init regions file bed_streamer, TODO move this code to Starling_run
        std::unique_ptr<bed_streamer> assemble_regions;
        if (! regions_file.empty())
        {
          const std::string bam_region(get_starling_bam_region_string(opt,dopt));
          bed_streamer *bedstr = new bed_streamer(regions_file.c_str(),bam_region.c_str());
          while ( bedstr->next() ){
            const bed_record *br = bedstr->get_record_ptr();
            known_pos_range2 range(br->begin,br->end);
            std::vector<std::string> words;
            boost::split(words,br->line,boost::is_any_of("\t"));
            std::vector<std::string> contigs(words.begin()+3,words.end());
            this->rt.addRegion(range,contigs);
          }
          delete bedstr;
        }
        else
        {
            //add in dummy dev regions
            known_pos_range2 range(239692924,239695935);
            this->rt.addRegion(range,std::vector<std::string>());
            known_pos_range2 range2(239692945,239695950);
            this->rt.addRegion(range2,std::vector<std::string>());
        }

    }

    bool keep_extending(int st, int end){
        return (this->rt.isPayloadInRegion(st) && this->rt.isPayloadInRegion(end));
    }
    bool do_assemble(int st, int end)
    {
        return (rt.isPayloadInRegion(st));
    }

    RegionPayloadTracker< std::vector< std::string > >  rt;

private:
    int assembleCount, assembleContigLength;          // count of regions to assemble, cummulative length of assembled regions
    std::string regions_file;
    const RegionTracker& _nocompress_regions;
    /// given an assembler with a region buffered, predict if it should it be assembled
    const starling_base_options& opt;
    const starling_deriv_options& dopt;
};

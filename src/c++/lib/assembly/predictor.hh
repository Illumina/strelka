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


namespace ASSEMBLY_TRIGGER
{

enum index_t
{
    bedTrack,
    highVarDensity,
	indelConflicts,
	NONE,
    SIZE
};

inline
const char*
get_label(const unsigned idx)
{
    switch (idx)
    {
    case bedTrack:
        return "BedTrack";
    case highVarDensity:
        return "highDensity";
    case indelConflicts:
        return "indelConflict";
    default:
        assert(false && "Unknown trigger");
        return nullptr;
    }
}
}

struct predictor
{
    predictor(
              const starling_base_options& init_opt,
              const starling_deriv_options& init_dopt,
              const RegionTracker& init_nocompress_regions)
      : regions_file(init_opt.assembly_regions_filename), opt(init_opt), dopt(init_dopt),
        _nocompress_regions(init_nocompress_regions),
        lengthToBreakRegion(20),
        varCountCutoff(3)

    {

        beginNotFullyProcessed = -1;
        beginPossibleAssemblyRegion = -1;
        numVarsInPossibleAssemblyRegion = 0;
        lastVarPos = -1;
        lastPos = -1;

        //init regions file bed_streamer, TODO move this code to Starling_run
        // JUNK: std::unique_ptr<bed_streamer> assemble_regions;

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
//        else
//        {
            //add in dummy dev regions
//            known_pos_range2 range(239691265,239691280);
//            this->rt.addRegion(range,std::vector<std::string>());
//            known_pos_range2 range2(239691282,239691285);
//            this->rt.addRegion(range2,std::vector<std::string>());
//        }

        // log_os << this->rt.regionCount() << std::endl;

    }

    void add_site( int pos, bool isVar )
    {

      if (beginNotFullyProcessed == -1)
      {
          beginNotFullyProcessed = pos;
      }

      if (isVar)
      {
          if (beginPossibleAssemblyRegion == -1)
          {
              beginPossibleAssemblyRegion = pos;
          }
          ++numVarsInPossibleAssemblyRegion;
          if(lastVarPos < pos)
          {
              lastVarPos = pos;
          }
      }
      if (lastPos < pos) // test this in case an indel already pushed us past here
      {
          lastPos = pos;
      }

#if 0
      if ( isVar )
      std::cerr << " site addition: "
                << beginNotFullyProcessed << " : "
                <<  beginPossibleAssemblyRegion << " : "
                << numVarsInPossibleAssemblyRegion << " : "
                << lastVarPos << " : "
                << lastPos << " : "
                <<  isVar << "\n";
#endif
    }



    void add_indel( int begpos, int endpos )
    {

      if (beginNotFullyProcessed == -1)
      {
          beginNotFullyProcessed = begpos;
      }

      if (beginPossibleAssemblyRegion == -1)
      {
          beginPossibleAssemblyRegion = begpos;
      }
      ++numVarsInPossibleAssemblyRegion;
      lastVarPos = endpos;

      lastPos = endpos;

#if 0
      std::cerr << " indel addition: "
                << beginNotFullyProcessed << " : "
                <<  beginPossibleAssemblyRegion << " : "
                << numVarsInPossibleAssemblyRegion << " : "
                << lastVarPos << " : "
                << lastPos << "\n";
#endif
    }

    bool keep_extending(){
        return ( lastPos-lastVarPos < lengthToBreakRegion);
    }

    known_pos_range2 do_assemble_and_update()
    {
        known_pos_range2 asmRange(-1,-1);
        if (numVarsInPossibleAssemblyRegion >= varCountCutoff)
        {


#if 0
      std::cerr << beginNotFullyProcessed << " : "
                <<  beginPossibleAssemblyRegion << " : "
                << numVarsInPossibleAssemblyRegion << " : "
                << lastVarPos << " : "
                << lastPos << "\n";
#endif

            asmRange.set_begin_pos(beginPossibleAssemblyRegion);
            asmRange.set_end_pos( lastVarPos+1 );
            assemblyReason = ASSEMBLY_TRIGGER::highVarDensity;
        }
        beginNotFullyProcessed = lastVarPos+1;
        numVarsInPossibleAssemblyRegion = 0;
        beginPossibleAssemblyRegion = -1;
        lastVarPos = -1;
        return asmRange;
    }

    RegionPayloadTracker< std::vector< std::string > >  rt;

    // predictor state
    int beginNotFullyProcessed;              // positions before here have already been fully processed: assembled or not assembled, but are no longer pending
    int beginPossibleAssemblyRegion;         // bases from here on might be involved in assembly
    int numVarsInPossibleAssemblyRegion;
    int lastVarPos;                          // last position that would encourage assembly
    int lastPos;                             // last position that has been incorporated into predictor so far

    // predictor parameters
    const int lengthToBreakRegion;           // if we see this many homref bases in a row, we know there is no assembly involving that interval
    const int varCountCutoff;                // if we see this many separate variations called without achieving lengthToBreakRegion, we should assemble

    ASSEMBLY_TRIGGER::index_t assemblyReason;
private:
    int assembleCount, assembleContigLength;          // count of regions to assemble, cummulative length of assembled regions
    std::string regions_file;
    const RegionTracker& _nocompress_regions;
    /// given an assembler with a region buffered, predict if it should it be assembled
    const starling_base_options& opt;
    const starling_deriv_options& dopt;
};

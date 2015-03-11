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
#include "blt_util/RegionTracker.hh"

//#define DEBUG_predictor


#ifdef DEBUG_predictor
#include "blt_util/log.hh"
#endif

struct predictor
{
    predictor(const RegionTracker& assembly_regions)
    	: _assembly_regions(assembly_regions)
    {

        //add in dummy dev regions
        known_pos_range2 range(239691265,239691280);
        this->rt.addRegion(range);
        known_pos_range2 range2(239691282,239691285);
        this->rt.addRegion(range2);
//        log_os << this->rt.regionCount() << std::endl;
    }
    bool keep_extending(int st, int end){
    	return (this->rt.isInRegion(st) && this->rt.isInRegion(end));
    }
    bool do_assemble(int st, int end)
    {
        return (this->rt.isInRegion(st));
    }

private:
    int assembleCount, assembleContigLength;          // count of regions to assemble, cummulative length of assembled regions
    std::string regions_file;
    RegionTracker rt;
    const RegionTracker& _assembly_regions;
    /// given an assembler with a region buffer, predict if it should be assembled
};

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
 *  Author: Morten Kallberg / Peter Krusche
 */

#pragma once

#include "applications/starling/site_info_stream.hh"
#include "applications/starling/gvcf_locus_info.hh"
#include "starling_common/starling_base_shared.hh"
#include "starling_common/starling_read_buffer.hh"
#include "starling_common/starling_pos_processor_util.hh"
#include "applications/starling/starling_shared.hh"

#include <list>

/**
 * @brief Predictor interface class
 * 
 * Defines that basic interface all predictors have to implement.
 * 
 * Implementations are hidden in predictor.cpp
 * 
 * The basic class 
 */
class predictor_interface
{
public:
    /** add sites / update positions */
//    virtual void add_site(site_info const & si) { }
//    virtual void add_indel(indel_info const & ii) { }

    /** keep extending the current block of variants
        (this tells assembler to not forward things on to 
         the gVCF aggregator)
     */
    virtual bool keep_extending() 
    {
        return false;
    }

    /**
     * @brief Return ranges to assemble
     * 
     * Last range will be [-1, -1)
     * 
     * Returned string can contain information in prediction reasons.
     */
    std::pair<known_pos_range2, std::string> next_range()
    {        
        return std::make_pair(known_pos_range2(-1, -1), std::string());
    }
};

enum class PredictorType { bedfile, density };

/**
 * @brief Predictor Factory
 * 
 * @param init_opt Starling command line options  
 * @param init_dopt Starling command line options 
 * @return new predictor interface (free using delete)
 * 
 * This might get more arguments as more predictors get implemented
 * (ideally, we'd do this using perfect forwarding, but then we can't
 *  link externally).
 * 
 */
extern std::unique_ptr<predictor_interface> make_predictor(
    PredictorType _type, 
    const starling_base_options& init_opt, 
    const starling_deriv_options& init_dopt);

/** split a stream using predictions */
class predictor_stream_splitter : public site_info_stream
{
public:
    ~predictor_stream_splitter()
    {
        // rather than flush put in assert() to check if buffers empty
        // 
        flush();
    }

    /** implement site_info_stream */
    bool add_site(site_info& si);
    bool add_indel(const indel_info& ii);
    void flush();

    /** add predictor 
     *  Maximally 64 predictors are supported as consumers are fed
     *  a mask 
     */
    uint64_t add_predictor(std::unique_ptr<predictor_interface> pred);

    /** add consumer 
     *  
     *  Predictor is chosen by matching either the predictor masks.
     *  
     *  predictor_mask_any_of gives a mask of predictors which can "say yes"
     *  to trigger this consumer.
     *  
     *  predictor_mask_none_of requires that none of the given predictors to "say yes"
     *  
     *  Variants that don't match any consumer get discarded. To add a default
     *  consumer use register_consumer()
     */
    void add_consumer(std::shared_ptr<site_info_stream> consumer,
                      uint64_t predictor_mask_any_of,
                      uint64_t predictor_mask_none_of = 0);

private:
    // list of predictors
    std::vector< std::unique_ptr<predictor_interface> > predictors; 

    // list of returned predictor segments with info (same size as predictors)
    typedef std::vector< std::pair<known_pos_range2, std::string> > seglist;    
    std::vector<seglist> seglists;

    // consumer list
    typedef struct _stream_output {
        uint64_t any; uint64_t none; std::shared_ptr<site_info_stream> consumer;
    } stream_output;

    std::list<stream_output> consumers; 
};

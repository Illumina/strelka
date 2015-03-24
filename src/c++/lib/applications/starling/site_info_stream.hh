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
 * site_info_stream.hh
 *
 *  Created on: Mar 10, 2015
 *      Author: Morten Kallberg
 *
 *  Abstract class specifying the site_info stream interface for
 *  linking together consumers of site_info / indel_infos
 */

#include "gvcf_locus_info.hh"

#include <deque>
#include <memory>

#ifndef C___SITE_INFO_STREAM_HH_
#define C___SITE_INFO_STREAM_HH_

class site_info_stream {
public:
    site_info_stream(){}

    virtual ~site_info_stream()
    {
        // TODO this is not super-tidy, it should be handled via 
        // flush() in derived class constructors.
        // Derived site_info_stream classes should clear the buffer in 
        // their destructor to make sure this doesn't have unwanted 
        // side-effects.
        notify_consumer_up_to();
    }

    /** interface methods */
    virtual bool add_site(site_info& si) = 0;
    virtual bool add_indel(const indel_info& ii) = 0;

    /** (PK) this should go away!
     *  ideally remove from interface, move up to caller code
     */
    virtual bool add_indel(const pos_t pos,
                   const indel_key ik,
                   const starling_diploid_indel_core& dindel,
                   const starling_indel_report_info& iri,
                   const starling_indel_sample_report_info& isri)
    {
        //TODO: is this inefficient enough to avoid the extra copy?
        indel_info ii;
        ii.init(pos,ik,dindel,iri,isri);
        return add_indel(ii);
    }

    /** notify and empty buffer */
    virtual void flush()  = 0;

    // TODO this needs to return a sensible value base on what is expected in starling_pos_processor around line 59
    virtual bool is_phasing_block(){ return true; }

    /** default consumer */
    void register_consumer(std::shared_ptr<site_info_stream> consumer)
    {
        this->_consumer = consumer;
    }

 protected:

    /** notify a single consumer up to a position, optionally clear buffers up to there */
    void notify_one_consumer(std::shared_ptr<site_info_stream> consumer, int stopPos, bool clear=false);

    /** notify default consumer, clear buffer -- stop pos < 0 means empty buffer */
    void notify_consumer_up_to(int stopPos = -1);

    /** forget buffer contents without notifying */
    void clear_buffer();

    /** forget partial buffer contents without notifying */
    void clear_site_buffer_up_to(std::deque<site_info>::iterator & endIt);
    void clear_site_buffer_to_pos(int stopPos);
    void clear_indel_buffer_up_to(std::deque<indel_info>::iterator & endIt);
    void clear_indel_buffer_to_pos(int stopPos);

protected:
    std::shared_ptr<site_info_stream> _consumer;
    std::deque<site_info> _site_buffer;
    std::deque<indel_info> _indel_buffer;
};

#endif /* C___SITE_INFO_STREAM_HH_ */

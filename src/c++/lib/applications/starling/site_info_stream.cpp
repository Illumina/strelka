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

/**
 * \brief Implementation of site_info_stream 
 *
 * \file site_info_stream.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "site_info_stream.hh"

void
site_info_stream::
notify_one_consumer(std::shared_ptr<site_info_stream> consumer, int stopPos, bool clear)
{
    std::deque<site_info>::iterator sit = _site_buffer.begin();
    std::deque<indel_info>::iterator iit = _indel_buffer.begin();
    while( (sit != _site_buffer.end()  && sit->pos < stopPos) 
        || (iit != _indel_buffer.end() && iit->pos < stopPos))
    {
        if (sit != _site_buffer.end())
        {
            if(iit != _indel_buffer.end())
            {
                if ((*sit).pos >= (*iit).pos)
                {
                    consumer->add_indel(*iit);
                    ++iit;
                    continue;
                }
            }
            consumer->add_site(*sit);
            ++sit;
            continue;
        }
        // iit != _indel_buffer.end() thanks to loop invariant
        consumer->add_indel(*iit);
        ++iit;
    }

    if(clear)
    {
        clear_site_buffer_up_to(sit);
        clear_indel_buffer_up_to(iit);
    }
}

void 
site_info_stream::
notify_consumer_up_to(int stopPos)
{
    if (!this->_consumer)
    {
        return;
    }

    if(stopPos < 0)
    {
        int lastSitePos = -1;
        if(_site_buffer.size()>0)
        {
            lastSitePos = _site_buffer.back().pos;
        }
        int lastIndelPos = -1;
        if(_indel_buffer.size()>0)
        {
            lastIndelPos = _indel_buffer.back().pos;
        }
        stopPos = std::max(lastSitePos,lastIndelPos) + 1;
    }
    notify_one_consumer(_consumer, stopPos, true);
}

void 
site_info_stream::
clear_buffer(){
    this->_indel_buffer.clear();
    this->_site_buffer.clear();
}


void 
site_info_stream::
clear_site_buffer_up_to(std::deque<site_info>::iterator & endIt){
    _site_buffer.erase(_site_buffer.begin(),endIt);
}

void 
site_info_stream::
clear_site_buffer_to_pos(int stopPos){
    std::deque<site_info>::iterator sit;
    for (sit =  this->_site_buffer.begin();sit < this->_site_buffer.end();++sit)
    {
        if ( sit->pos >= stopPos)
        {
            break;
        }
    }
    clear_site_buffer_up_to(sit);
}

void 
site_info_stream::
clear_indel_buffer_up_to(std::deque<indel_info>::iterator & endIt){
    _indel_buffer.erase(_indel_buffer.begin(),endIt);
}

void 
site_info_stream::
clear_indel_buffer_to_pos(int stopPos){
    std::deque<indel_info>::iterator iit;
    for (iit =  this->_indel_buffer.begin();iit < this->_indel_buffer.end();++iit)
    {
        if ( iit->pos >= stopPos)
        {
            break;
        }
    }
    clear_indel_buffer_up_to(iit);
}

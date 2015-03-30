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
 *  predictor.cpp
 *
 *  Created on: Sep 10, 2013
 *  Author: Morten Kallberg / Peter Krusche
 */

#include "predictor.hh"

#include "predictor_bed.hh"
#include "predictor_density.hh"

template<typename T, typename... Args>
static std::unique_ptr<T> make_unique(Args&&... args)
{
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

/**
 * @brief Predictor Factory
 * 
 * @param init_opt Starling command line options  
 * @param init_dopt Starling command line options 
 * @return new predictor interface (free using delete)
 */
std::unique_ptr<predictor_interface> 
make_predictor(
    PredictorType _type,
    const starling_base_options& init_opt, const starling_deriv_options& init_dopt)
{
    switch(_type)
    {
        case PredictorType::bedfile:
            return make_unique<predictor_bed>(init_opt.assembly_regions_filename.c_str(), 
                                              get_starling_bam_region_string(init_opt, 
                                                                             init_dopt).c_str());
        case PredictorType::density:
            return make_unique<predictor_density>();
        default:
            assert(false && "Unknown trigger");
            break;
    }
}

bool 
predictor_stream_splitter::
add_site(site_info& si)
{
    bool success = true;
    bool still_collecting = false;
    for(auto & predictor : predictors)
    {
        predictor->add_site(si);
        still_collecting = still_collecting || predictor->keep_extending();
    }

    if(!still_collecting)
    {
        flush();
    }

    return success;
}

bool 
predictor_stream_splitter::
add_indel(const indel_info& ii)
{
    bool success = true;
    bool still_collecting = false;
    for(auto & predictor : predictors)
    {
        predictor->add_indel(ii);
        still_collecting = still_collecting || predictor->keep_extending();
    }

    if(!still_collecting)
    {
        flush();
    }

    return success;
}

void 
predictor_stream_splitter::
flush()
{
    int i = 0;
    size_t total_segs = 0;
    int min_pos = -1;

    for(auto & predictor : predictors)
    {
        std::pair<known_pos_range2, std::string> p = predictor->next_range();
        while(p.first.begin_pos() >= 0)
        {
            seglists[i].push_back(p);
            if(min_pos < 0)
            {
                min_pos = p.first.begin_pos();
            }
            else
            {
                min_pos = std::min(p.first.begin_pos(), min_pos);
            }
            ++total_segs;
            p = predictor->next_range();
        }
        ++i;
    }
    // elementary intervals
    int max_pos = 1;

    while(total_segs && max_pos > 0)
    {
        // drop sites not covered by any predictor
        if(_consumer)
        {
            notify_one_consumer(_consumer, min_pos + 1, true);
            _consumer->flush();
        }
        else
        {
            std::deque<site_info>::iterator sit = _site_buffer.begin();
            while( sit != _site_buffer.end() && sit->pos < min_pos)
            {
                ++sit;
            }
            std::deque<indel_info>::iterator iit = _indel_buffer.begin();
            while( iit != _indel_buffer.end() && iit->pos < min_pos)
            {
                ++iit;
            }
            clear_site_buffer_up_to(sit);
            clear_indel_buffer_up_to(iit);
        }

        for(auto & s : seglists)
        {
            while(!s.empty())
            {
                if(s.front().first.end_pos() < min_pos + 1)
                {
                    s.erase(s.begin());
                    --total_segs;
                }
            }
        }

        uint64_t current = 0;
        uint64_t cm = 1;
        max_pos = -1;

        for(auto const & s : seglists)
        {
            if(!s.empty() && s.front().first.begin_pos() == min_pos)
            {
                current |= cm;

                if(max_pos < 0)
                {
                    max_pos = s.front().first.end_pos();
                }
                else
                {
                    max_pos = std::min(s.front().first.end_pos(), max_pos);
                }
            }
            cm <<= 1;
        }

        // TODO determine if no consumer notified, then send to default
        // rewrite such that only single consumer. flush on consumer switch.
        for (auto & output : consumers)
        {
            if((output.any & current) && ((output.none & current) == 0))
            {
                notify_one_consumer(output.consumer, max_pos+1, false);
                // TODO collect all consumers which received variants, 
                // flush ones that haven't received at the end.
                output.consumer->flush();
            }
        }
        min_pos = max_pos;
    }

    if(_consumer)
    {
        int lastSitePos = -1;
        if(_site_buffer.size() > 0)
        {
            lastSitePos = _site_buffer.back().pos;
        }
        int lastIndelPos = -1;
        if(_indel_buffer.size() > 0)
        {
            lastIndelPos = _indel_buffer.back().pos;
        }

        notify_one_consumer(_consumer, std::max(lastSitePos, lastIndelPos) + 1, true);
        _consumer->flush();
    }

    _site_buffer.clear();
    _indel_buffer.clear();
}

/** add predictor 
 *  Maximally 64 predictors are supported as consumers are fed
 *  a mask 
 */
uint64_t 
predictor_stream_splitter::
add_predictor(std::unique_ptr<predictor_interface> pred)
{
    assert (predictors.size() < 64 && "Maximally 64 predictors can be added");
    uint64_t mask = 1 << (int8_t)predictors.size();
    predictors.push_back(std::move(pred));
    seglists.resize(predictors.size());
    return mask;
}

/** add consumer 
 *  
 *  Predictor is chosen by matching either the predictor masks.
 *  
 *  predictor_mask_any_of gives a mask of predictors which can "say yes"
 *  to trigger this consumer.
 *  
 *  predictor_mask_all_of requires all given predictors to "say yes"
 */
void 
predictor_stream_splitter::
add_consumer(
    std::shared_ptr<site_info_stream> consumer, 
    uint64_t predictor_mask_any_of,
    uint64_t predictor_mask_none_of)
{
    consumers.push_back(stream_output{predictor_mask_any_of, predictor_mask_none_of, consumer});
}

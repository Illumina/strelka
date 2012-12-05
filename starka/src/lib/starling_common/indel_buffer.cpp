// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2012 Illumina, Inc.
//
// This software is covered by the "Illumina Genome Analyzer Software
// License Agreement" and the "Illumina Source Code License Agreement",
// and certain third party copyright/licenses, and any user of this
// source file is bound by the terms therein (see accompanying files
// Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
// Illumina_Source_Code_License_Agreement.pdf and third party
// copyright/license notices).
//
//

/// \file
///
/// \author Chris Saunders
///

#include "starling_common/indel_buffer.hh"

#include <cassert>

#include <iostream>


//#define IB_DEBUG



std::pair<indel_buffer::iterator,indel_buffer::iterator>
indel_buffer::
pos_range_iter(const pos_t begin_pos,
               const pos_t end_pos) {
    const indel_key end_range_key(end_pos);
    const iterator end(_idata.lower_bound(end_range_key));
    const indel_key begin_range_key(begin_pos-static_cast<pos_t>(_max_indel_size));
    iterator begin(_idata.lower_bound(begin_range_key));
    for(;begin!=end;++begin) {
        if(begin->first.right_pos() >= begin_pos) break;
    }
    return std::make_pair(begin,end);
}



// The goal is to return all indels with a left or right breakpoint in the
// range [begin_pos,end_pos]. Returning indels in addition to this set is
// acceptable.
//
// The indels_keys: "end_range_key" and "begin_range_key" take
// advantage of the indel NONE type, which sorts ahead of all other
// types at the same position.
//
std::pair<indel_buffer::const_iterator,indel_buffer::const_iterator>
indel_buffer::
pos_range_iter(const pos_t begin_pos,
               const pos_t end_pos) const {
    const indel_key end_range_key(end_pos);
    const const_iterator end(_idata.lower_bound(end_range_key));
    const indel_key begin_range_key(begin_pos-static_cast<pos_t>(_max_indel_size));
    const_iterator begin(_idata.lower_bound(begin_range_key));
    for(;begin!=end;++begin) {
        if(begin->first.right_pos() >= begin_pos) break;
    }
    return std::make_pair(begin,end);
}



bool
indel_buffer::
insert_indel(const indel_observation& obs,
             const bool is_shared,
             bool& is_repeat_obs) {

    assert(obs.key.type != INDEL::NONE);
    idata_t::iterator i(_idata.find(obs.key));
    if(i == _idata.end()){
        indel_data id;
        id.add_observation(obs.data,is_shared,is_repeat_obs);
        _idata.insert(std::make_pair(obs.key,id));
        return true;
    }

    indel_data& id(get_indel_data(i));
    id.add_observation(obs.data,is_shared,is_repeat_obs);
    return false;
}



void
indel_buffer::
clear_pos(const pos_t pos) {
    const iterator i_begin(pos_iter(pos));
    const iterator i_end(pos_iter(pos+1));
    _idata.erase(i_begin,i_end);
}


static
void
dump_range(indel_buffer::const_iterator i,
           const indel_buffer::const_iterator i_end,
           std::ostream& os)
{
    for(;i!=i_end;++i){
        os << i->first << get_indel_data(i);
    }
}


void
indel_buffer::
dump_pos(const pos_t pos,
         std::ostream& os) const
{
    dump_range(pos_iter(pos),pos_iter(pos+1),os);
}



void
indel_buffer::
dump(std::ostream& os) const
{
    os << "INDEL_BUFFER DUMP ON\n";
    dump_range(_idata.begin(),_idata.end(),os);
    os << "INDEL_BUFFER DUMP OFF\n";
}

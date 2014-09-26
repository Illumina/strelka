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

///
/// \author Chris Saunders
///

#include "starling_common/pos_basecall_buffer.hh"

#include <iostream>



// debug dumpers:
void
pos_basecall_buffer::
dump(std::ostream& /*os*/) const
{
#if 0
    for (const auto& val : _pdata)
    {
        os << "pc_buff pos: " << val.first << "\n";
    }
#endif
}

void
pos_basecall_buffer::
update_ranksums(
    char refpos,
    const pos_t pos,
    const base_call& bc,
    const uint8_t adjustedMapq,
    const unsigned cycle)
{
    const bool is_reference(refpos==id_to_base(bc.base_id));

    auto& posdata(_pdata.getRef(pos));
    posdata.baseq_ranksum.add_observation(is_reference,static_cast<unsigned>(bc.get_qscore()));
    posdata.mq_ranksum.add_observation(is_reference,static_cast<unsigned>(adjustedMapq));
    posdata.read_pos_ranksum.add_observation(is_reference,cycle);
}


void
pos_basecall_buffer::
update_read_pos_ranksum(
    char refchar,
    const pos_t pos,
    const base_call& bc,
    const unsigned read_pos)
{
    const bool is_reference(refchar==id_to_base(bc.base_id));

    auto& posdata(_pdata.getRef(pos));
    posdata.read_pos_ranksum.add_observation(is_reference,read_pos);
}

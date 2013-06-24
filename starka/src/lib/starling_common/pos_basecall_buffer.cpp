// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/downloads/sequencing/licenses/>.
//

/// \file
///
/// \author Chris Saunders
///

#include "starling_common/pos_basecall_buffer.hh"

#include <iostream>



// debug dumpers:
void
pos_basecall_buffer::
dump(std::ostream& os) const {

    pciter i(_pdata.begin()), i_end(_pdata.end());
    for(; i!=i_end; ++i) {
        os << "pc_buff pos: " << i->first << "\n";
    }
}

void
pos_basecall_buffer::update_ranksums(char refpos, const pos_t pos,const base_call& bc, const uint8_t mapq, const int cycle){
	const bool is_reference(refpos==id_to_base(bc.base_id));

	_pdata[pos].baseq_ranksum.add_observation(is_reference,static_cast<int>(bc.get_qscore()));
	_pdata[pos].mq_ranksum.add_observation(is_reference,static_cast<int>(mapq));
	_pdata[pos].read_pos_ranksum.add_observation(is_reference,cycle);
}

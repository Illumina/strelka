//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2018 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
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
updateGermlineScoringMetrics(
    const char refchar,
    const pos_t pos,
    const uint8_t call_id,
    const uint8_t qscore,
    const uint8_t mapq,
    const unsigned cycle,
    const unsigned distanceFromReadEdge,
    const bool is_submapped)
{
    const bool is_reference(refchar==id_to_base(call_id));

    auto& posdata(_pdata.getRef(pos));
    posdata.mq_ranksum.add_observation(is_reference,static_cast<unsigned>(mapq));
    if (! is_submapped)
    {
        posdata.baseq_ranksum.add_observation(is_reference,static_cast<unsigned>(qscore));
        posdata.readPositionRankSum.add_observation(is_reference,cycle);
        if (not is_reference)
        {
            // maxDistance may help the mean edge distance better generalize over different read lengths:
            static const unsigned maxDistanceFromEdge(20);
            posdata.distanceFromReadEdge.addObservation(std::min(maxDistanceFromEdge, distanceFromReadEdge));
        }
    }
}


void
pos_basecall_buffer::
update_read_pos_ranksum(
    char refchar,
    const pos_t pos,
    const uint8_t call_id,
    const unsigned read_pos)
{
    const bool is_reference(refchar==id_to_base(call_id));

    auto& posdata(_pdata.getRef(pos));
    posdata.readPositionRankSum.add_observation(is_reference,read_pos);
}

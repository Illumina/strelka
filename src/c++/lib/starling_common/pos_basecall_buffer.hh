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

#pragma once

#include "blt_common/snp_pos_info.hh"
#include "blt_util/blt_types.hh"
#include "blt_util/RangeMap.hh"

#include <iosfwd>
#include <cmath>
#include <string>


struct EmptyPosSet
{
    EmptyPosSet()
    {
        for (unsigned i(0); i<BASE_ID::SIZE; ++i)
        {
            pis[i].set_ref_base(id_to_base(i));
        }
    }
    std::array<snp_pos_info,BASE_ID::SIZE> pis;
};



struct pos_basecall_buffer
{
    pos_basecall_buffer(
        const reference_contig_segment& ref)
        : _ref(ref), _pdata(ref)
    {}

    void
    clear()
    {
        _pdata.clear();
    }

    void
    insert_pos_submap_count(const pos_t pos)
    {
        _pdata.getRef(pos).submappedReadCount++;
    }

    void
    insert_pos_spandel_count(const pos_t pos)
    {
        _pdata.getRef(pos).spanningDeletionReadCount++;
    }

    // update mapQ sum for RMSMappingQuality calculation
    void
    insert_mapq_count(
        const pos_t pos,
        const uint8_t mapq)
    {
        _pdata.getRef(pos).mapqTracker.add(mapq);
    }

    void
    insert_alt_read_pos(
        const pos_t pos,
        const uint8_t call_id,
        const uint16_t readPos,
        const uint16_t readLength)
    {
        snp_pos_info& posdata(_pdata.getRef(pos));
        if (posdata.get_ref_base() == id_to_base(call_id)) return;

        posdata.altAlleleReadPositionInfo.push_back({readPos,readLength});
    }

    /// add associated basecall data to pileup at position pos which
    /// can be used to compute variant scoring metrics
    ///
    void
    updateGermlineScoringMetrics(
        const char refpos,
        const pos_t pos,
        const uint8_t call_id,
        const uint8_t qscore,
        const uint8_t mapq,
        const unsigned cycle,
        const unsigned distanceFromReadEdge,
        const bool is_submapped);

    void
    update_read_pos_ranksum(
        char refchar,
        const pos_t pos,
        const uint8_t call_id,
        const unsigned read_pos);

    void
    insert_pos_basecall(const pos_t pos,
                        const bool is_tier1,
                        const base_call& bc)
    {
        if (is_tier1)
        {
            _pdata.getRef(pos).calls.push_back(bc);
        }
        else
        {
            _pdata.getRef(pos).tier2_calls.push_back(bc);
        }
    }

    void
    decrementSpanningIndelPloidy(const pos_t pos)
    {
        _pdata.getRef(pos).spanningIndelPloidyModification -= 1;
    }

    const snp_pos_info&
    get_pos(const pos_t pos) const
    {
        static const EmptyPosSet emptyPosSet;
        if (! _pdata.isKeyPresent(pos))
        {
            return emptyPosSet.pis[base_to_id(_ref.get_base(pos))];
        }
        return _pdata.getConstRef(pos);
    }

    void
    clear_pos(const pos_t pos)
    {
        if (_pdata.isKeyPresent(pos)) _pdata.erase(pos);
    }

    /// clear all data up to and including pos
    void
    clear_to_pos(const pos_t pos)
    {
        _pdata.eraseTo(pos);
    }

    bool
    empty() const
    {
        return _pdata.empty();
    }

    void
    dump(std::ostream& os) const;

private:
    typedef RangeMap<pos_t,snp_pos_info,ClearT<snp_pos_info>> pdata_t;

    // inherit so that we can intercept the getRef calls:
    struct PosData : public pdata_t
    {
        PosData(const reference_contig_segment& ref_init) : ref(ref_init) {}

        snp_pos_info&
        getRef(
            const pos_t& pos)
        {
            snp_pos_info& pi(pdata_t::getRef(pos));
            if (! pi.is_ref_set()) pi.set_ref_base(ref.get_base(pos));
            return pi;
        }

        const reference_contig_segment& ref;
    };

    const reference_contig_segment& _ref;
    PosData _pdata;
};


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
        for (unsigned i(0);i<BASE_ID::SIZE;++i)
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
    insert_pos_submap_count(const pos_t pos)
    {
        _pdata.getRef(pos).n_submapped++;
    }

    void
    insert_pos_spandel_count(const pos_t pos)
    {
        _pdata.getRef(pos).n_spandel++;
    }

    // update mapQ sum for MQ calculation
    void
    insert_mapq_count(
        const pos_t pos,
        const uint8_t mapq,
        const uint8_t adjustedMapq)
    {
        snp_pos_info& posdata(_pdata.getRef(pos));
        posdata.n_mapq++;
        posdata.cumm_mapq += (adjustedMapq*adjustedMapq);
        if (mapq==0) posdata.n_mapq0++;
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

        posdata.altReadPos.push_back({readPos,readLength});
    }

    // add single base meta-data to rank-sum pile-up data-structures
    void
    update_ranksums(
        char refpos,
        const pos_t pos,
        const uint8_t call_id,
        const uint8_t qscore,
        const uint8_t adjustedMapq,
        const unsigned cycle,
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
    insert_hap_cand(const pos_t pos,
                    const bool is_tier1,
                    const bam_seq_base& read_seq,
                    const uint8_t* qual,
                    const unsigned offset)
    {
        // TODO write this for multi-tier:
        assert(is_tier1);
        _pdata.getRef(pos).hap_set.emplace_back(read_seq,qual,offset);
    }

    const snp_pos_info&
    get_pos(const pos_t pos) const
    {
        static const EmptyPosSet empty;
        if (! _pdata.isKeyPresent(pos))
        {
            return empty.pis[base_to_id(_ref.get_base(pos))];
        }
        return _pdata.getConstRef(pos);
    }

    void
    clear_pos(const pos_t pos)
    {
        if (_pdata.isKeyPresent(pos)) _pdata.erase(pos);
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


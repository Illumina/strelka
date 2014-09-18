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
#include "blt_util/rangeMap.hh"

#include <iosfwd>
#include <cmath>
#include <string>


struct pos_basecall_buffer
{
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

    // add single base meta-data to rank-sum pile-up data-structures
    void
    update_ranksums(
        char refpos,
        const pos_t pos,
        const base_call& bc,
        const uint8_t adjustedMapq,
        const unsigned cycle);

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

    // returns NULL for empty pos
    //
    snp_pos_info*
    get_pos(const pos_t pos)
    {
        if (! _pdata.isKeyPresent(pos)) return nullptr;
        return &_pdata.getRef(pos);
    }

    const snp_pos_info*
    get_pos(const pos_t pos) const
    {
        if (! _pdata.isKeyPresent(pos)) return nullptr;
        return &_pdata.getConstRef(pos);
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

#if 0
    iterator pos_iter(const pos_t pos)
    {
        return _pdata.lower_bound(pos);
    }
    const_iterator pos_iter(const pos_t pos) const
    {
        return _pdata.lower_bound(pos);
    }

    // debug dumpers:
    void
    dump_pos(const pos_t pos, std::ostream& os) const;
#endif

    void
    dump(std::ostream& os) const;

private:
    typedef rangeMap<pos_t,snp_pos_info,ClearT<snp_pos_info>> pdata_t;

    pdata_t _pdata;
};


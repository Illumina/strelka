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

#pragma once

#include "blt_common/snp_pos_info.hh"
#include "blt_util/blt_types.hh"

#include <iosfwd>
#include <map>
#include <cmath>
#include <string>


struct pos_basecall_buffer {

    void
    insert_pos_submap_count(const pos_t pos) {

        _pdata[pos].n_submapped++;
    }

    void
    insert_pos_spandel_count(const pos_t pos) {

        _pdata[pos].n_spandel++;
    }

    // update mapQ sum for MQ calculation
    void
    insert_mapq_count(const pos_t pos, const uint8_t mapq) {

        _pdata[pos].n_mapq++;
        _pdata[pos].cumm_mapq += (mapq*mapq);
        //we calculate the RMS, so store squared mapq
    }

    // add single base meta-data to rank-sum pile-up data-structures
    void
    update_ranksums(char refpos, const pos_t pos,const base_call& bc, const uint8_t mapq, const int cycle);

    void
    insert_pos_basecall(const pos_t pos,
                        const bool is_tier1,
                        const base_call& bc) {
        if(is_tier1) {
            _pdata[pos].calls.push_back(bc);
        } else {
            _pdata[pos].tier2_calls.push_back(bc);
        }
    }

    void
    insert_hap_cand(const pos_t pos,
                    const bool is_tier1,
                    const bam_seq_base& read_seq,
                    const uint8_t* qual,
                    const unsigned offset) {

        // TODO write this for multi-tier:
        assert(is_tier1);
        _pdata[pos].hap_set.push_back(hap_cand(read_seq,qual,offset));
    }

    // returns NULL for empty pos
    //
    snp_pos_info*
    get_pos(const pos_t pos) {
        const piter i(_pdata.find(pos));
        if(i==_pdata.end()) return NULL;
        return &(i->second);
    }

    const snp_pos_info*
    get_pos(const pos_t pos) const {
        const pciter i(_pdata.find(pos));
        if(i==_pdata.end()) return NULL;
        return &(i->second);
    }

    void
    clear_pos(const pos_t pos) {
        const piter i(_pdata.find(pos));
        if(i!=_pdata.end()) _pdata.erase(i);
    }

    bool
    empty() const { return _pdata.empty(); }


#if 0
    iterator pos_iter(const pos_t pos) {
        return _pdata.lower_bound(pos);
    }
    const_iterator pos_iter(const pos_t pos) const {
        return _pdata.lower_bound(pos);
    }

    // debug dumpers:
    void
    dump_pos(const pos_t pos, std::ostream& os) const;
#endif

    void
    dump(std::ostream& os) const;

private:
    typedef std::map<pos_t,snp_pos_info> pdata_t;
    typedef pdata_t::iterator piter;
    typedef pdata_t::const_iterator pciter;

    pdata_t _pdata;
};


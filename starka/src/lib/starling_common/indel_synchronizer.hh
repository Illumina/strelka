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

#ifndef __INDEL_SYNCHRONIZER_HH
#define __INDEL_SYNCHRONIZER_HH

#include "blt_util/id_map.hh"
#include "starling_common/depth_buffer.hh"
#include "starling_common/indel_buffer.hh"
#include "starling_common/starling_shared.hh"


struct indel_synchronizer;


struct indel_sync_data {

    void
    register_sample(indel_buffer& ib,
                    const depth_buffer& db,
                    const starling_sample_options& sample_opt,
                    const sample_id_t sample_no);

private:
    friend struct indel_synchronizer;

    struct indel_sample_data {

        indel_sample_data()
            : ibp(0)
            , dbp(0)
            , sample_optp(0)
        {}

        indel_sample_data(indel_buffer& ib,
                          const depth_buffer& db,
                          const starling_sample_options& sample_opt)
            : ibp(&ib)
            , dbp(&db)
            , sample_optp(&sample_opt)
        {}

        indel_buffer* ibp;
        const depth_buffer* dbp;
        const starling_sample_options* sample_optp;
    };

    typedef id_map<sample_id_t,indel_sample_data> idata_t;

    idata_t _idata;
};



// this structure helps to sync the indel information from multiple
// samples, currently used for tumor/normal indel-calling.
//
// There is one indel synchronizer associated with each sample
// (referred to as the primary sample below). The synchronizer defines
// the primary sample's synchronization policy with any other sample.
//
struct indel_synchronizer {

    // ctor for simple single-sample operation:
    //
    indel_synchronizer(indel_buffer& ib,
                       const depth_buffer& db,
                       const starling_sample_options& init_sample_opt)
        : _sample_no(0)
        , _sample_order(0)
    {
        _isd.register_sample(ib,db,init_sample_opt,_sample_no);
    }

    // ctor for multi-sample synced cases:
    //
    // sample_no is the sample that is 'primary'
    // for this syncronizer.
    //
    indel_synchronizer(const indel_sync_data& isd,
                       const sample_id_t sample_no)
        : _isd(isd)
        , _sample_no(sample_no)
        , _sample_order(_isd._idata.get_id(sample_no)) {}


    indel_buffer&
    ibuff() { return ibuff(_sample_order); }

    const indel_buffer&
    ibuff() const { return ibuff(_sample_order); }


    // returns true if this indel is novel to the buffer
    //
    // indel is fully inserted into the primary sample buffer, but
    // only the key is inserted into other sample buffers.
    bool
    insert_indel(const indel_observation& obs);

    // is an indel treated as a candidate for genotype calling and
    // realignment or as a "private" (ie. noise) indel?
    //
    bool
    is_candidate_indel(const starling_options& opt,
                       const indel_key& ik,
                       const indel_data& id) const {

        if(not id.status.is_candidate_indel_cached) {
            is_candidate_indel_int(opt,ik,id);
        }
        return id.status.is_candidate_indel;
    }

    // this version is less efficient than if you have indel_data
    // beforehand, but provided for convenience:
    //
    bool
    is_candidate_indel(const starling_options& opt,
                       const indel_key& ik) const {
        const indel_data* id_ptr(ibuff().get_indel_data_ptr(ik));
        assert(NULL != id_ptr);
        return is_candidate_indel(opt,ik,*id_ptr);
    }

    // used for debug output:
    sample_id_t
    get_sample_id() const {
        return _sample_no;
    }

private:

    void
    is_candidate_indel_int(const starling_options& opt,
                           const indel_key& ik,
                           const indel_data& id) const;

    indel_buffer&
    ibuff(const unsigned s) { return *(idata().get_value(s).ibp); }
    const indel_buffer&
    ibuff(const unsigned s) const { return *(idata().get_value(s).ibp); }

    const depth_buffer&
    ebuff(const unsigned s) const { return *(idata().get_value(s).dbp); }

    const starling_sample_options&
    sample_opt(const unsigned s) const { return *(idata().get_value(s).sample_optp); }

    typedef indel_sync_data::idata_t idata_t;

    idata_t&
    idata() { return _isd._idata; }

    const idata_t&
    idata() const { return _isd._idata; }

    indel_sync_data _isd;

    // this is the "external" id of the primary sample, it can be
    // considered as a map key
    //
    const sample_id_t _sample_no;

    // this is the "internal" sequential id of the primary sample, it
    // can be considered an array index
    //
    const unsigned _sample_order;
};


#endif

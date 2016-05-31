// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2016 Illumina, Inc.
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

#include "min_count_binom_gte_cache.hh"

#include "blt_util/depth_buffer.hh"
#include "blt_util/id_map.hh"
#include "starling_common/indel_buffer.hh"
#include "starling_common/starling_base_shared.hh"

#include <vector>


struct indel_sample_data
{
    indel_sample_data(
        indel_buffer& ib,
        const depth_buffer& db,
        const depth_buffer& db2,
        const starling_sample_options& sample_opt,
        const double init_max_depth)
        : ibp(&ib)
        , dbp(&db)
        , dbp2(&db2)
        , sample_optp(&sample_opt)
        , max_depth(init_max_depth)
    {}

    indel_buffer* ibp;
    const depth_buffer* dbp;
    const depth_buffer* dbp2;
    const starling_sample_options* sample_optp;
    double max_depth;
};


/// helps to sync the indel information from multiple
/// samples, currently used for tumor/normal indel-calling.
///
/// There is one indel synchronizer associated with each sample
/// (referred to as the primary sample below). The synchronizer defines
/// the primary sample's synchronization policy with any other sample.
///
struct indel_synchronizer
{
    /// instantiate for multi-sample synced cases:
    ///
    /// \param[in] sample_no is the sample that is 'primary' for this synchronizer.
    ///
    indel_synchronizer(
        const starling_base_options& opt,
        const starling_base_deriv_options& dopt,
        const reference_contig_segment& ref)
        : _opt(opt)
        , _dopt(dopt)
        , _ref(ref)
        , _countCache(dopt.countCache)
    {}

    unsigned
    register_sample(
        indel_buffer& ib,
        const depth_buffer& db,
        const depth_buffer& db2,
        const starling_sample_options& sample_opt,
        const double max_depth);

    void
    finalizeSamples()
    {
        assert(! _idata.empty());
        _isFinalized = true;
    }

    indel_buffer&
    ibuff(
        const unsigned sampleId)
    {
        return *idata(sampleId).ibp;
    }

    const indel_buffer&
    ibuff(
        const unsigned sampleId) const
    {
        return *idata(sampleId).ibp;
    }

    /// \returns true if this indel is novel to the buffer
    ///
    /// indel is fully inserted into the primary sample buffer, but
    /// only the key is inserted into other sample buffers.
    bool
    insert_indel(
        const unsigned sampleId,
        const indel_observation& obs);

    /// is an indel treated as a candidate for genotype calling and
    /// realignment or as a "private" (ie. noise) indel?
    ///
    bool
    is_candidate_indel(
        const unsigned sampleId,
        const indel_key& ik,
        const indel_data& id) const
    {
        if (! id.status.is_candidate_indel_cached)
        {
            is_candidate_indel_impl(sampleId, ik, id);
        }
        return id.status.is_candidate_indel;
    }

    // this version is less efficient than if you have indel_data
    // beforehand, but provided for convenience:
    //
    bool
    is_candidate_indel(
        const unsigned sampleId,
        const indel_key& ik) const
    {
        const indel_data* id_ptr(ibuff(sampleId).get_indel_data_ptr(ik));
        if (nullptr == id_ptr) find_data_exception(ik);
        return is_candidate_indel(sampleId, ik, *id_ptr);
    }

private:

    /// test whether the indel should be promoted to
    /// a candidate
    bool
    is_candidate_indel_impl_test_signal_noise(
        const indel_key& ik,
        const indel_data& id,
        const indel_data* idsp[],
        const unsigned isds) const;

    /// much weaker version of the above -- used for indel
    /// discovery protocols outside of the variant caller
    bool
    is_candidate_indel_impl_test_weak_signal(
        const indel_data* idsp[],
        const unsigned isds) const;

    bool
    is_candidate_indel_impl_test(
        const indel_key& ik,
        const indel_data& id,
        const indel_data* idsp[],
        const unsigned sampleCount) const;

    void
    is_candidate_indel_impl(
        const unsigned sampleId,
        const indel_key& ik,
        const indel_data& id) const;

    /// return object which provides estimated depth of tier1 reads
    const depth_buffer&
    ebuff(const unsigned sampleId) const
    {
        return *(idata(sampleId).dbp);
    }

    /// return object which provides estimated depth of tier2 reads
    const depth_buffer&
    ebuff2(const unsigned sampleId) const
    {
        return *(idata(sampleId).dbp2);
    }

    const starling_sample_options&
    sample_opt(
        const unsigned sampleId) const
    {
        return *(idata(sampleId).sample_optp);
    }

    unsigned
    getSampleCount() const
    {
        return _idata.size();
    }

    typedef std::vector<indel_sample_data> idata_t;

    indel_sample_data&
    idata(
        const unsigned sampleId)
    {
        assert(_isFinalized);
        assert(sampleId<_idata.size());

        return _idata[sampleId];
    }

    const indel_sample_data&
    idata(
        const unsigned sampleId) const
    {
        assert(_isFinalized);
        assert(sampleId<_idata.size());

        return _idata[sampleId];
    }

    void
    find_data_exception(const indel_key& ik) const;


    const starling_base_options& _opt;
    const starling_base_deriv_options& _dopt;
    const reference_contig_segment& _ref;
    const min_count_binom_gte_cache& _countCache;

    bool _isFinalized = false;
    idata_t _idata;
};

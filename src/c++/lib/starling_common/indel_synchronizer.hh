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


#include "blt_util/depth_buffer.hh"
#include "starling_common/indel.hh"
#include "starling_common/min_count_binom_gte_cache.hh"
#include "starling_common/starling_base_shared.hh"

#include <vector>


/// coordinate indel information across multiple samples
///
struct indel_synchronizer
{
    indel_synchronizer(
        const starling_base_options& opt,
        const starling_base_deriv_options& dopt,
        const reference_contig_segment& ref)
        : _opt(opt)
        , _dopt(dopt)
        , _ref(ref)
        , _countCache(dopt.countCache)
    {}

    /// prior to executing any other functions, each sample must be registered:
    unsigned
    register_sample(
        const depth_buffer &db,
        const depth_buffer &db2,
        const starling_sample_options &sample_opt,
        const double max_depth);

    /// registration must be followed with a finalization step:
    void
    finalizeSamples()
    {
        assert(! _indelSampleData.empty());
        _isFinalized = true;
    }

    typedef IndelData indel_buffer_value_t;
    typedef std::map<indel_key,indel_buffer_value_t> indel_buffer_t;
    typedef indel_buffer_t::iterator iterator;
    typedef indel_buffer_t::const_iterator const_iterator;

    /// \returns true if this indel is novel to the buffer
    ///
    bool
    insert_indel(
        const unsigned sampleId,
        const indel_observation& obs);

    /// position iterators based on left-most indel position:
    iterator
    pos_iter(const pos_t pos)
    {
        return _indelBuffer.lower_bound(indel_key(pos));
    }

    const_iterator
    pos_iter(const pos_t pos) const
    {
        return _indelBuffer.lower_bound(indel_key(pos));
    }

    /// position iterators which return (at least) all indels with a
    /// left or right breakpoint in the range.
    ///
    /// Note:
    /// 1) indels which encompass the range are not returned
    /// 2) some non-intersecting indels may be returned in the
    ///    iteration range
    ///
    std::pair<iterator,iterator>
    pos_range_iter(const pos_t begin_pos, const pos_t end_pos);

    std::pair<const_iterator,const_iterator>
    pos_range_iter(const pos_t begin_pos, const pos_t end_pos) const;

    /// return nullptr if no indel found:
    indel_buffer_value_t*
    get_indel_data_ptr(const indel_key& ik)
    {
        const iterator i(_indelBuffer.find(ik));
        return ((i==_indelBuffer.end()) ? nullptr : &(i->second) );
    }

    const indel_buffer_value_t*
    get_indel_data_ptr(const indel_key& ik) const
    {
        const const_iterator i(_indelBuffer.find(ik));
        return ((i==_indelBuffer.end()) ? nullptr : &(i->second) );
    }

    /// is an indel treated as a candidate for genotype calling and
    /// realignment or as a "private" (ie. noise) indel?
    ///
    bool
    is_candidate_indel(
        const indel_key& ik,
        const IndelData& id) const
    {
        if (! id.status.is_candidate_indel_cached)
        {
            is_candidate_indel_impl(ik, id);
        }
        return id.status.is_candidate_indel;
    }

    /// this version is less efficient than if you have indel_data
    /// beforehand, but provided for convenience:
    ///
    bool
    is_candidate_indel(
        const indel_key& ik) const
    {
        const IndelData* id_ptr(get_indel_data_ptr(ik));
        if (nullptr == id_ptr) find_data_exception(ik);
        return is_candidate_indel(ik, *id_ptr);
    }

    void
    clear_pos(const pos_t pos);

    bool
    empty() const
    {
        return _indelBuffer.empty();
    }

    // debug dumpers:
    void
    dump_pos(const pos_t pos, std::ostream& os) const;

    void
    dump(std::ostream& os) const;

private:

    /// helper struct for indel_syncronizer
    struct indel_sample_data
    {
        indel_sample_data(
            const depth_buffer& db,
            const depth_buffer& db2,
            const starling_sample_options& sample_opt,
            const double init_max_depth)
            : dbp(&db)
            , dbp2(&db2)
            , sample_optp(&sample_opt)
            , max_depth(init_max_depth)
        {}

        const depth_buffer* dbp;
        const depth_buffer* dbp2;
        const starling_sample_options* sample_optp;
        double max_depth;
    };

    /// test whether the indel should be promoted to
    /// a candidate
    bool
    is_candidate_indel_impl_test_signal_noise(
        const indel_key& ik,
        const IndelData& id) const;

    /// much weaker version of the above -- used for indel
    /// discovery protocols outside of the variant caller
    bool
    is_candidate_indel_impl_test_weak_signal(
        const IndelData& id) const;

    bool
    is_candidate_indel_impl_test(
        const indel_key& ik,
        const IndelData& id) const;

    void
    is_candidate_indel_impl(
        const indel_key& ik,
        const IndelData& id) const;

    /// return object which provides estimated depth of tier1 reads
    const depth_buffer&
    ebuff(const unsigned sampleId) const
    {
        return *(getIndelSampleData(sampleId).dbp);
    }

    /// return object which provides estimated depth of tier2 reads
    const depth_buffer&
    ebuff2(const unsigned sampleId) const
    {
        return *(getIndelSampleData(sampleId).dbp2);
    }

    const starling_sample_options&
    sample_opt(
        const unsigned sampleId) const
    {
        return *(getIndelSampleData(sampleId).sample_optp);
    }

    unsigned
    getSampleCount() const
    {
        return _indelSampleData.size();
    }

    typedef std::vector<indel_sample_data> indel_sample_data_t;

    indel_sample_data&
    getIndelSampleData(
        const unsigned sampleId)
    {
        assert(_isFinalized);
        assert(sampleId<_indelSampleData.size());
        return _indelSampleData[sampleId];
    }

    const indel_sample_data&
    getIndelSampleData(
        const unsigned sampleId) const
    {
        assert(_isFinalized);
        assert(sampleId<_indelSampleData.size());
        return _indelSampleData[sampleId];
    }

    void
    find_data_exception(const indel_key& ik) const;

/////////////// data
    const starling_base_options& _opt;
    const starling_base_deriv_options& _dopt;
    const reference_contig_segment& _ref;
    const min_count_binom_gte_cache& _countCache;

    bool _isFinalized = false;
    indel_sample_data_t _indelSampleData;
    indel_buffer_t _indelBuffer;
};



// These functions assume valid iterators:
//
// note these are just stopgaps so that client code can
// voluntarily abstract out the storage mechanism of indel_data,
// ideally we would not give raw iterators to client code, see
// top-level TODO
inline
indel_synchronizer::indel_buffer_value_t&
get_indel_data(const indel_synchronizer::iterator i)
{
    return (i->second);
}

inline
const indel_synchronizer::indel_buffer_value_t&
get_indel_data(const indel_synchronizer::const_iterator i)
{
    return (i->second);
}


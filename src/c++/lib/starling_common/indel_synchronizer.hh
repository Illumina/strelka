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

# include <vector>

#include "blt_util/id_map.hh"
#include "starling_common/depth_buffer.hh"
#include "starling_common/indel_buffer.hh"
#include "starling_common/starling_base_shared.hh"


struct indel_synchronizer;


struct indel_sync_data
{
    void
    register_sample(indel_buffer& ib,
                    const depth_buffer& db,
                    const depth_buffer& db2,
                    const starling_sample_options& sample_opt,
                    const double max_depth,
                    const sample_id_t sample_no);

private:
    friend struct indel_synchronizer;

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

    typedef id_map<sample_id_t,indel_sample_data> idata_t;

    idata_t _idata;
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
    /// instantiate for simple single-sample operation:
    ///
    /// \param[in] max_candidate_depth - max depth (in this sample) for indel candidates, any filtration will be applied to all samples. A negative value disables the filter.
    ///
    /// \max_candidate_depth - max depth (in this sample) for indel candidates, any filtration will be applied to all samples. A negative value disables the fi
    ///
    indel_synchronizer(
        const starling_base_options& opt,
        const reference_contig_segment& ref,
        const double max_candidate_depth,
        indel_buffer& ib,
        const depth_buffer& db,
        const depth_buffer& db2,
        const starling_sample_options& init_sample_opt)
        : _opt(opt)
        , _ref(ref)
        , _sample_no(0)
        , _sample_order(0)
    {
        _isd.register_sample(ib,db,db2,init_sample_opt, max_candidate_depth,_sample_no);
    }

    /// instantiate for multi-sample synced cases:
    ///
    /// \param[in] sample_no is the sample that is 'primary' for this synchronizer.
    ///
    indel_synchronizer(
        const starling_base_options& opt,
        const reference_contig_segment& ref,
        const indel_sync_data& isd,
        const sample_id_t sample_no)
        : _opt(opt)
        , _ref(ref)
        , _isd(isd)
        , _sample_no(sample_no)
        , _sample_order(_isd._idata.get_id(sample_no)) {}

    indel_buffer&
    ibuff()
    {
        return ibuff(_sample_order);
    }

    const indel_buffer&
    ibuff() const
    {
        return ibuff(_sample_order);
    }

    /// \returns true if this indel is novel to the buffer
    ///
    /// indel is fully inserted into the primary sample buffer, but
    /// only the key is inserted into other sample buffers.
    bool
    insert_indel(const indel_observation& obs);

    /// is an indel treated as a candidate for genotype calling and
    /// realignment or as a "private" (ie. noise) indel?
    ///
    bool
    is_candidate_indel(
        const indel_key& ik,
        const indel_data& id) const
    {
        if (! id.status.is_candidate_indel_cached)
        {
            is_candidate_indel_impl(ik,id);
        }
        return id.status.is_candidate_indel;
    }

    // this version is less efficient than if you have indel_data
    // beforehand, but provided for convenience:
    //
    bool
    is_candidate_indel(
        const indel_key& ik) const
    {
        const indel_data* id_ptr(ibuff().get_indel_data_ptr(ik));
        if (nullptr == id_ptr) find_data_exception(ik);
        return is_candidate_indel(ik,*id_ptr);
    }

    // used for debug output:
    sample_id_t
    get_sample_id() const
    {
        return _sample_no;
    }

private:

    bool
    is_candidate_indel_impl_test(
        const indel_key& ik,
        const indel_data& id,
        const indel_data* idsp[],
        const unsigned isds) const;

    void
    is_candidate_indel_impl(
        const indel_key& ik,
        const indel_data& id) const;

    indel_buffer&
    ibuff(const unsigned s)
    {
        return *(idata().get_value(s).ibp);
    }

    const indel_buffer&
    ibuff(const unsigned s) const
    {
        return *(idata().get_value(s).ibp);
    }

    /// return object which provides estimated depth of tier1 reads
    const depth_buffer&
    ebuff(const unsigned sample) const
    {
        return *(idata().get_value(sample).dbp);
    }

    /// return object which provides estimated depth of tier2 reads
    const depth_buffer&
    ebuff2(const unsigned sample) const
    {
        return *(idata().get_value(sample).dbp2);
    }

    const starling_sample_options&
    sample_opt(const unsigned s) const
    {
        return *(idata().get_value(s).sample_optp);
    }

    typedef indel_sync_data::idata_t idata_t;

    idata_t&
    idata()
    {
        return _isd._idata;
    }

    const idata_t&
    idata() const
    {
        return _isd._idata;
    }

    void
    find_data_exception(const indel_key& ik) const;


    const starling_base_options& _opt;
    const reference_contig_segment& _ref;

    indel_sync_data _isd;

    // Rather than calculate p-values explicitly in the hpol one case
    // from a binomial distribution, or use a lookup table for total cov
    // x indel size, we can leverage the fact the n * p ~ n * p * (1 - p)
    // for our current indel error rates and use a Poisson approximation
    // to the binomial distribution.   This will not work when p is large
    // relative to n (i.e. high error rates--in a typical 30x dataset, error
    // rates would have to be around 10%, we're currently observing error rates
    // around 5e-6 - 2e-5.  This lookup table provides the minimum coverage
    // given an error rate p and coverage n when for a p-value of 1e-9.
    // We could calculate this dynamically in a small piece of code while we
    // do refCoverage, etc. for an arbitrary p-value, but this should do for now.

    /* Generated from the following code in R
    *   require(plyr)

    *   min_success_enum <- data.frame(min_success = qpois(1e-9, lambda = seq(1e-6, 3, 1e-6), 
    *                                                      lower.tail = FALSE),
    *                                  np = seq(1e-6, 3, 1e-6))
    *   min_success_enum$min_success[min_success_enum$min_success < 2] <- 2

    *   min_success_ranges <- ddply(min_success_enum, .(min_success), summarize, 
    *                               min_success = unique(min_success),
    *                               min_np = min(np),
    *                               max_np = max(np))
    *   min_success_ranges$min_np <- min_success_ranges$min_np - 1e-6

    *   cat(apply(min_success_ranges, 1, function(x) sprintf("{%6d, %.6f, %.6f},",
    *                                             x["min_success"],
    *                                             x["min_np"], x["max_np"])),
    *       sep = "\n")

    */
    // min success min_np    max_np
    static constexpr const double _min_hpol_one_indel_cov_by_np[17][3] = {
        {     2, 0.000000, 0.001817},
        {     3, 0.001817, 0.012477},
        {     4, 0.012477, 0.041576},
        {     5, 0.041576, 0.095978},
        {     6, 0.095978, 0.179024},
        {     7, 0.179024, 0.291566},
        {     8, 0.291566, 0.433032},
        {     9, 0.433032, 0.602135},
        {    10, 0.602135, 0.797280},
        {    11, 0.797280, 1.016790},
        {    12, 1.016790, 1.259019},
        {    13, 1.259019, 1.522411},
        {    14, 1.522411, 1.805519},
        {    15, 1.805519, 2.107018},
        {    16, 2.107018, 2.425694},
        {    17, 2.425694, 2.760443},
        {    18, 2.760443, 3.000000}
    };

    static constexpr const unsigned _min_hpol_one_indel_cov_rows = 17;
    static_assert(starling_base_options::tumor_min_hpol_pval == 1e-9,
        "Lookup table is invalid for minimum candidate p-values != 1e-9.  Please refer to src/c++/lib/starling_common/indel_synchronizer.hh:264, recalculate and replace lookup table, then recompile");
        
    unsigned get_min_candidate_cov(double np) const
    {
        // this should be faster than a binary interval search since
        // we almost always expect n * p to be on the low side
        // (i.e. np < 0.01)
        for(unsigned i = 0; i < _min_hpol_one_indel_cov_rows; ++i)
        {
            if(np <= _min_hpol_one_indel_cov_by_np[i][2])
            {
                return (unsigned )_min_hpol_one_indel_cov_by_np[i][0];
            }
        }
        return 0U;
    }

    // this is the "external" id of the primary sample, it can be
    // considered as a map key
    //
    const sample_id_t _sample_no;

    // this is the "internal" sequential id of the primary sample, it
    // can be considered an array index
    //
    const unsigned _sample_order;
};

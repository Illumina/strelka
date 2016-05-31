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
/// \author Mitch Bekritsky
///

#include "starling_common/indel_synchronizer.hh"

#include "blt_util/binomial_test.hh"
#include "blt_util/blt_exception.hh"
#include "blt_util/log.hh"
#include "calibration/IndelErrorModel.hh"

#include <cstdlib>

#include <iostream>
#include <sstream>



unsigned
indel_synchronizer::
register_sample(
    indel_buffer& ib,
    const depth_buffer& db,
    const depth_buffer& db2,
    const starling_sample_options& sample_opt,
    const double max_depth)
{
    assert(! _isFinalized);
    const unsigned sampleIndex(_idata.size());
    _idata.emplace_back(ib,db,db2,sample_opt,max_depth);
    return sampleIndex;
}



bool
indel_synchronizer::
insert_indel(
    const unsigned sampleId,
    const indel_observation& obs)
{
    // first insert indel into this sample:
    bool is_synced_sample(false);
    bool is_repeat_obs(false);

    const bool is_novel(ibuff(sampleId).insert_indel(obs,is_synced_sample,is_repeat_obs));

    // then insert indel into synchronized samples:
    is_synced_sample=true;
    const unsigned sampleCount(getSampleCount());
    for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
    {
        if (sampleId) continue;
        ibuff(sampleIndex).insert_indel(obs,is_synced_sample,is_repeat_obs);
    }
    return is_novel;
}



bool
indel_synchronizer::
is_candidate_indel_impl_test_signal_noise(
    const indel_key& ik,
    const indel_data& id,
    const indel_data* idsp[],
    const unsigned isds) const
{
    //
    // Step 1: find the error rate for this indel and context:
    //
    double indelToRefErrorProb(0.);
    {
        starling_indel_report_info iri;
        get_starling_indel_report_info(ik,id,_ref,iri);

        // refToIndelErrorProb does not factor in to this calculation, so put this in a temp value
        double refToIndelErrorProb(0.);
        _dopt.getIndelErrorModel().getIndelErrorRate(iri,refToIndelErrorProb,indelToRefErrorProb);
    }

    //
    // Step 2: determine if the observed counts are sig wrt the error rate for at least one sample:
    //
    for (unsigned i(0); i<isds; ++i)
    {
        // for each sample, get the number of tier 1 reads supporting the indel
        // and the total number of tier 1 reads at this locus
        const unsigned n_indel_reads = idsp[i]->all_read_ids.size();
        unsigned n_total_reads = ebuff(i).val(ik.pos-1);

        // total reads and indel reads are measured in different ways here, so the nonsensical
        // result of indel_reads>total_reads is possible. The total is fudged below to appear sane
        // before going into the count test:
        n_total_reads = std::max(n_total_reads,n_indel_reads);

        if (n_total_reads == 0) continue;

        // this min indel support coverage is applied regardless of error model:
        static const unsigned min_candidate_cov_floor(2);
        if (n_indel_reads < min_candidate_cov_floor) continue;

        // test to see if the observed indel coverage has a binomial exact test
        // p-value above the rejection threshold. If this does not occur for the
        // counts observed in any sample, the indel cannot become a candidate
        if (_countCache.isRejectNull(n_total_reads, indelToRefErrorProb, n_indel_reads))
        {
            return true;
        }
    }
    return false;
}



bool
indel_synchronizer::
is_candidate_indel_impl_test_weak_signal(
    const indel_data* idsp[],
    const unsigned isds) const
{
    for (unsigned i(0); i<isds; ++i)
    {
        // for each sample, get the number of tier 1 reads supporting the indel
        const unsigned n_indel_reads = idsp[i]->all_read_ids.size();

        static const unsigned min_candidate_cov_floor(1);
        if (n_indel_reads < min_candidate_cov_floor) continue;

        return true;
    }
    return false;
}



bool
indel_synchronizer::
is_candidate_indel_impl_test(
    const indel_key& ik,
    const indel_data& id,
    const indel_data* idsp[],
    const unsigned sampleCount) const
{
    // check whether the candidate has been externally specified:
    for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
    {
        if (idsp[sampleIndex]->is_external_candidate) return true;
    }

    if (_opt.is_candidate_indel_signal_test)
    {
        if (! is_candidate_indel_impl_test_signal_noise(ik,id,idsp,sampleCount)) return false;
    }
    else
    {
        if (! is_candidate_indel_impl_test_weak_signal(idsp,sampleCount)) return false;
    }

    /////////////////////////////////////////
    // test against short open-ended segments:
    //
    // call get_insert_size() instead of asking for the insert seq
    // so as to not finalize any incomplete insertions:
    if (ik.is_breakpoint() &&
        (_opt.min_candidate_indel_open_length > id.get_insert_size()))
    {
        return false;
    }

    /////////////////////////////////////////
    // test against max_depth:
    //
    for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
    {
        const double max_depth(idata(sampleIndex).max_depth);
        if (max_depth <= 0.) continue;

        const unsigned estdepth(ebuff(sampleIndex).val(ik.pos-1));
        const unsigned estdepth2(ebuff2(sampleIndex).val(ik.pos-1));
        if ((estdepth+estdepth2) > max_depth) return false;
    }

    return true;
}



void
indel_synchronizer::
is_candidate_indel_impl(
    const unsigned sampleId,
    const indel_key& ik,
    const indel_data& id) const
{
    //////////////////////////////////////
    // lookup all indel_data objects:
    //
    const indel_data* idsp[MAX_SAMPLE];

    const unsigned sampleCount(getSampleCount());
    for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
    {
        if (sampleIndex==sampleId)
        {
            idsp[sampleIndex] = &id;
        }
        else
        {
            idsp[sampleIndex] = ibuff(sampleIndex).get_indel_data_ptr(ik);
            assert(nullptr != idsp[sampleIndex]);
        }
    }

    const bool is_candidate(is_candidate_indel_impl_test(ik,id,idsp,sampleCount));

    for (unsigned i(0); i<sampleCount; ++i)
    {
        idsp[i]->status.is_candidate_indel=is_candidate;
        idsp[i]->status.is_candidate_indel_cached=true;
    }
}



void
indel_synchronizer::
find_data_exception(const indel_key& ik) const
{
    std::ostringstream oss;
    oss << "ERROR: could not find indel_data for indel: " << ik;
    throw blt_exception(oss.str().c_str());
}


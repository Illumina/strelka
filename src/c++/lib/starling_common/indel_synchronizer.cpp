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
#include "starling_common/starling_indel_error_prob.hh"
#include "calibration/scoringmodels.hh"

#include <cstdlib>

#include <iostream>
#include <sstream>
#include <vector>



void
indel_sync_data::
register_sample(
    indel_buffer& ib,
    const depth_buffer& db,
    const depth_buffer& db2,
    const starling_sample_options& sample_opt,
    const double max_depth,
    const sample_id_t sample_no)
{
    if (_idata.test_key(sample_no))
    {
        log_os << "ERROR: sample_no " << sample_no << " repeated in indel sync registration\n";
        exit(EXIT_FAILURE);
    }
    _idata.insert(sample_no,indel_sample_data(ib,db,db2,sample_opt,max_depth));
}



bool
indel_synchronizer::
insert_indel(const indel_observation& obs)
{
    // first insert indel into this sample:
    bool is_synced_sample(false);
    bool is_repeat_obs(false);

    const bool is_novel(ibuff(_sample_order).insert_indel(obs,is_synced_sample,is_repeat_obs));

    // then insert indel into synchronized samples:
    is_synced_sample=true;
    const unsigned isds(idata().size());
    for (unsigned i(0); i<isds; ++i)
    {
        if (i == _sample_order) continue;
        ibuff(i).insert_indel(obs,is_synced_sample,is_repeat_obs);
    }
    return is_novel;
}


bool
indel_synchronizer::
is_candidate_indel_impl_test(
    const indel_key& ik,
    const indel_data& id,
    const indel_data* idsp[],
    const unsigned isds) const
{
    // check whether the candidate has been externally specified:
    for (unsigned i(0); i<isds; ++i)
    {
        if (idsp[i]->is_external_candidate) return true;
    }

    // have we already tested if the indel is much more likely to be noise than sample variation?
    bool is_indel_noise_checked(false);

    const IndelErrorModel indel_model = scoring_models::Instance().get_indel_model();

    // since we don't use indel size to calculate error rates at the moment in candidcacy
    // we can retrieve this value once
    static const indel_error_rates hpol_one_error_rates = indel_model.calc_abstract_prop(1, 1, 1, false);


    starling_indel_report_info iri;
    get_starling_indel_report_info(ik,id,_ref,iri);

    bool is_indel(ik.type == INDEL::INSERT || ik.type == INDEL::DELETE);

    if (is_indel)
    {
        bool is_pass_indel_noise_check(false);
        // HPOL ONE CASE
        // we're only testing against the reference error rate here
        // in hpol one cases, we need to establish a few things:
        // 1. reference tract length is 0 or 1
        // 2. indel error rate
        // 3. total coverage
        if (iri.ref_repeat_count <= 1)
        {
            is_indel_noise_checked = true;

            for (unsigned i(0); i<isds; ++i)
            {
                const unsigned n_indel_reads = idsp[i]->all_read_ids.size();
                unsigned n_total_reads = ebuff(i).val(ik.pos-1);
                n_total_reads = std::max(n_total_reads,n_indel_reads);

                if (n_total_reads == 0)
                {
                    continue;
                }

                // test to see if the observed indel coverage exceeds the minimum value
                // for the rejection threshold.

                // this min indel support coverage is applied regardless of error model:
                static const unsigned min_candidate_cov_floor(2);
                if (n_indel_reads < min_candidate_cov_floor)
                {
                    continue;
                }

                // this min indel support coverage is based on error model:
                double error_rate(0.);
                if     (ik.type == INDEL::INSERT)
                {
                    error_rate = hpol_one_error_rates.insert_rate;
                }
                else
                {
                    error_rate = hpol_one_error_rates.delete_rate;
                }
                is_pass_indel_noise_check=_countCache.isRejectNull(n_total_reads, error_rate, n_indel_reads);

                // if any sample passes the noise check, the indel is a candidate
                if (is_pass_indel_noise_check)
                {
                    break;
                }
            }

            if (! is_pass_indel_noise_check) return false;
        }

        // NON-HPOL ONE CASE
        // check broader STR criteria for non-hpol one
        if (!is_indel_noise_checked)
        {
            is_indel_noise_checked=true;
            bool use_length_dependence=false;

            double ref_error_prob(0.);
            double indel_error_prob(0.);
            // indel_error_prob does not factor in to this calculation, but is
            // required by get_indel_error_prob

            // get expected per-read error rate for this STR
            indel_model.calc_prop(_opt,iri,indel_error_prob,ref_error_prob,use_length_dependence);

            // need to undo the reference-error bias term that is used to reduce overcalls
            ref_error_prob /= _opt.indel_ref_error_factor;

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

                // test to see if the observed indel coverage has a binomial exact test
                // p-value below the rejection threshold.  If it does not, it cannot be
                // a candidate indel
                is_pass_indel_noise_check=(is_reject_binomial_gte_n_success_exact(_opt.tumor_min_hpol_pval, ref_error_prob,
                                                                                  n_indel_reads, n_total_reads));

                // if any sample passes the homopolymer noise check (i.e. indel has high enough
                // coverage that it may not be noise), the indel is a candidate
                if (is_pass_indel_noise_check)
                {
                    break;
                }
            }

            if (! is_pass_indel_noise_check) return false;
        }
    }

    assert(!is_indel || is_indel_noise_checked);

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
    for (unsigned i(0); i<isds; ++i)
    {
        const double max_depth(idata().get_value(i).max_depth);
        if (max_depth <= 0.) continue;

        const unsigned estdepth(ebuff(i).val(ik.pos-1));
        const unsigned estdepth2(ebuff2(i).val(ik.pos-1));
        if ((estdepth+estdepth2) > max_depth) return false;
    }

    return true;
}



void
indel_synchronizer::
is_candidate_indel_impl(
    const indel_key& ik,
    const indel_data& id) const
{
    //////////////////////////////////////
    // lookup all indel_data objects:
    //
    const indel_data* idsp[MAX_SAMPLE];

    const unsigned isds(idata().size());
    for (unsigned i(0); i<isds; ++i)
    {
        if (i==_sample_order)
        {
            idsp[i] = &id;
        }
        else
        {
            idsp[i] = ibuff(i).get_indel_data_ptr(ik);
            assert(nullptr != idsp[i]);
        }
    }

    const bool is_candidate(is_candidate_indel_impl_test(ik,id,idsp,isds));

    for (unsigned i(0); i<isds; ++i)
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


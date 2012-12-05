// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2012 Illumina, Inc.
//
// This software is covered by the "Illumina Genome Analyzer Software
// License Agreement" and the "Illumina Source Code License Agreement",
// and certain third party copyright/licenses, and any user of this
// source file is bound by the terms therein (see accompanying files
// Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
// Illumina_Source_Code_License_Agreement.pdf and third party
// copyright/license notices).
//
//

/// \file
///
/// \author Chris Saunders
///

#include "blt_util/log.hh"
#include "starling_common/indel_synchronizer.hh"

#include <cstdlib>
#include <iostream>



void
indel_sync_data::
register_sample(indel_buffer& ib,
                const depth_buffer& db,
                const starling_sample_options& sample_opt,
                const sample_id_t sample_no) {

    if(_idata.test_key(sample_no)) {
        log_os << "ERROR: sample_no " << sample_no << " repeated in indel sync registration\n";
        exit(EXIT_FAILURE);
    }
    _idata.insert(sample_no,indel_sample_data(ib,db,sample_opt));
}



bool
indel_synchronizer::
insert_indel(const indel_observation& obs) {

    // first insert indel into this sample:
    bool is_synced_sample(false);
    bool is_repeat_obs(false);
    const bool is_novel(ibuff(_sample_order).insert_indel(obs,is_synced_sample,is_repeat_obs));

    // then insert indel into synchronized samples:
    is_synced_sample=true;
    const unsigned isds(idata().size());
    for(unsigned i(0);i<isds;++i) {
        if (i == _sample_order) continue;
        ibuff(i).insert_indel(obs,is_synced_sample,is_repeat_obs);
    }
    return is_novel;
}




void
indel_synchronizer::
is_candidate_indel_int(const starling_options& opt,
                       const indel_key& ik,
                       const indel_data& id) const {

    //////////////////////////////////////
    // lookup all indel_data objects:
    //
    const indel_data* idsp[MAX_SAMPLE];

    const unsigned isds(idata().size());
    for(unsigned i(0);i<isds;++i) {
        if(i==_sample_order) {
            idsp[i] = &id;
        } else {
            idsp[i] = ibuff(i).get_indel_data_ptr(ik);
            assert(NULL != idsp[i]);
        }
    }

    // pre-set result to false until candidacy is shown:
    for(unsigned i(0);i<isds;++i) {
        idsp[i]->status.is_candidate_indel=false;
        idsp[i]->status.is_candidate_indel_cached=true;
    }

    // check whether the candidate has been externally specified:
    bool is_external_candidate=false;
    for(unsigned i(0);i<isds;++i) {
        if(idsp[i]->is_external_candidate) {
            is_external_candidate=true;
            break;
        }
    }

    if(is_external_candidate) {
        for(unsigned i(0);i<isds;++i) {
            idsp[i]->status.is_candidate_indel=true;
        }
        return;
    }

    //////////////////////////////////////
    // test against min read count:
    //
    {
        bool is_min_count(false);

        int n_total_reads(0);
        for(unsigned i(0);i<isds;++i) {
            const int n_reads(idsp[i]->all_read_ids.size());

            // do the candidate reads exceed the (possibly lower than
            // default) sample specific threshold?:
            if(n_reads >= sample_opt(i).min_candidate_indel_reads) {
                is_min_count=true;
                break;
            }
            n_total_reads += n_reads;
        }

        // do reads from all samples exceed the default threshold?:
        if(n_total_reads >= opt.default_min_candidate_indel_reads) is_min_count=true;

        if(! is_min_count) return;
    }

    //////////////////////////////////////
    // test against min read frac:
    //
    {
        bool is_min_frac(false);

        double min_large_indel_frac(opt.min_candidate_indel_read_frac);
        const bool is_small_indel(static_cast<int>(std::max(ik.length,ik.swap_dlength)) <= opt.max_small_candidate_indel_size);

        // this value is used to get around type-mismatch error in
        // std::max() below
        static const unsigned one(1);

        for(unsigned i(0);i<isds;++i) {
            // note estdepth is based on genomic reads only, so
            // readfrac can be > 1:
            //
            const unsigned n_reads(idsp[i]->all_read_ids.size());
            const unsigned estdepth(std::max(one,ebuff(i).val(ik.pos-1)));
            const double readfrac(static_cast<double>(n_reads)/static_cast<double>(estdepth));

            double min_indel_frac(min_large_indel_frac);
            if(is_small_indel) {
                min_indel_frac=std::max(min_indel_frac,sample_opt(i).min_small_candidate_indel_read_frac);
            }

            // min_frac threshold only needs to pass in one sample to
            // be a candidate in all synchronized samples:
            if(readfrac >= min_indel_frac) {
                is_min_frac=true;
                break;
            }
        }

        if(! is_min_frac) return;
    }

    /////////////////////////////////////////
    // test against max_depth:
    //
    {
        const double max_depth(opt.max_candidate_indel_depth);
        for(unsigned i(0);i<isds;++i) {
            const unsigned estdepth(ebuff(i).val(ik.pos-1));
            if(estdepth > max_depth) return;
        }
    }

    /////////////////////////////////////////
    // test against short open-ended segments:
    //
    {
        if(ik.is_breakpoint() &&
           (opt.min_candidate_indel_open_length > id.get_insert_seq().size())) {
            return;
        }
    }

    // made it!
    for(unsigned i(0);i<isds;++i) {
        idsp[i]->status.is_candidate_indel=true;
    }
}

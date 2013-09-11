// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

/// \file
///
/// \author Chris Saunders
///

#pragma once

#include "blt_util/ranksum.hh"
#include "starling_common/indel_align_type.hh"
#include "starling_common/indel_key.hh"
#include "starling_common/starling_types.hh"

#include <cassert>
#include <iosfwd>
#include <map>
#include <string>
#include <set>
#include <vector>

//#define ID_DEBUG

#ifdef ID_DEBUG
#include "blt_util/log.hh"

#include <iostream>
#endif


/// represents the data associated with a single indel observation:
///
struct indel_observation_data {

    indel_observation_data()
        : is_noise(false)
        , is_external_candidate(false)
        , is_forced_output(false)
        , iat(INDEL_ALIGN_TYPE::GENOME_SUBMAP_READ)
        , id(0)
    {}

    bool is_noise;
    bool is_external_candidate;
    bool is_forced_output; // results of gt tests must be output even for very unlikely cases
    INDEL_ALIGN_TYPE::index_t iat;
    align_id_t id;
    std::string insert_seq;
};



// Holds the alignment scores created by each read aligning across an
// indel which express the relative probability of the read aligning
// to the indel or the reference (or elsewhere in the genome). Used
// for indel genotyping.
//
struct read_path_scores {

    typedef float score_t;

    read_path_scores(const score_t r=0,
                     const score_t i=0,
                     const uint16_t ns=0,
                     const uint16_t rlen=0,
                     const bool is_t1=true)
        : ref(r)
        , indel(i)
        , nsite(ns)
        , read_length(rlen)
        , is_tier1_read(is_t1)
    {}

    void
    insert_alt(const indel_key& ik,
               const score_t a) {
        const unsigned ais(static_cast<unsigned>(alt_indel.size()));
        if (ais < 2) {
            alt_indel.push_back(std::make_pair(ik,a));
        } else {
            unsigned min_index(ais);
            score_t min(a);
            for (unsigned i(0); i<ais; ++i) {
                if (alt_indel[i].second < min) {
                    min = alt_indel[i].second;
                    min_index = i;
                }
            }
            if (min_index<ais) {
                alt_indel[min_index] = std::make_pair(ik,a);
            }
        }
    }

    score_t ref;
    score_t indel;
    uint16_t nsite;

//    score_t alt;

    // store up to 2 highest scoring alternate indels
    typedef std::vector<std::pair<indel_key,score_t> > alt_indel_t;
    alt_indel_t alt_indel;

    // used to set expected het allele ratio:
    uint16_t read_length;

    // used to filter for/against tier2 data:
    bool is_tier1_read;
};



/// Accumulates evidence of the consensus insert sequence.
/// These classes make no attempt to find two insertions of the
/// same length with different insert sequences, but for these
/// cases this class will at least insure that we test and report
/// the most common insert sequence.
///
/// Note that for open-breakends we also report the longest insert
/// candidate, and consider the most prevalent as a secondary term.
///
struct insert_seq_manager {

    insert_seq_manager()
        : _is_consensus(false)
        , _obs_count(0)
    {}

    // get final consensus sequence:
    const std::string&
    get() {
        if (! _is_consensus) {
            _finalize();
        }
        return _consensus_seq;
    }

    // add insert sequence observation:
    void
    add_obs(const std::string& seq) {
        if (_is_consensus) {
            _exception("Attempting to add insert observation after finalizing");
        }

        // if we don't know what the most common insert will be after
        // this many samples are collected, then more cases will be
        // unlikely to help:
        static const unsigned max_obs_count(256);
        if (_obs_count>max_obs_count) return;

        obs_t::iterator i(_obs.find(seq));
        if (i == _obs.end()) {
            _obs[seq] = 1;
        } else {
            i->second += 1;
        }
        _obs_count += 1;
    }

private:
    void _exception(const char* msg) const;

    void _finalize();


    typedef std::map<std::string,unsigned> obs_t;
    bool _is_consensus;
    std::string _consensus_seq;
    unsigned _obs_count;
    obs_t _obs;
};



// represents the data from all observations associated with an indel
//
struct indel_data {

    indel_data(const indel_key& ik)
        : _ik(ik),
          is_external_candidate(false),
          is_forced_output(false)
    {}

    /// add an observation for this indel
    ///
    /// is_repeat_obs - has this read_id been observed before? this is both
    ///                 read and set by this method. Read ids are allowed to be
    ///                 repeated due to suggested alternate alignments from
    ///                 GROUPER
    void
    add_observation(const indel_observation_data& obs_data,
                    const bool is_shared,
                    bool& is_repeat_obs) {


#ifdef ID_DEBUG
        log_os << "KATTER: adding obs for indel: " << _ik;
        log_os << "KATTER: is_shared: " << is_shared << " is_repeat: " << is_repeat_obs << "\n";
        log_os << "KATTER: is_external: " << obs_data.is_external_candidate << " align_id: " << obs_data.id << "\n\n";
#endif

        if (! is_shared) {
            add_observation_core(obs_data,is_repeat_obs);
        }

        if (! obs_data.insert_seq.empty()) {
            if (! (is_shared && is_repeat_obs)) {
                _insert_seq.add_obs(obs_data.insert_seq);
            }
        }
    }




#if 0
    // add read and contig id evidence from another indel_data
    // structure -- this is not a copy ctor
    //
    void
    add_indel_data_evidence(const indel_data& id) {
        add_evidence(contig_ids,id.contig_ids);
        add_evidence(all_read_ids,id.all_read_ids);
        add_evidence(tier2_map_read_ids,id.tier2_map_read_ids);
        add_evidence(submap_read_ids,id.submap_read_ids);
        add_evidence(noise_read_ids,id.noise_read_ids);
        if (id.is_external_candidate) is_external_candidate=true;
    }
#endif


    const std::string&
    get_insert_seq() const {
#ifdef ID_DEBUG
        log_os << "KATTER: reporting insert seq for indel: " << _ik << "\n";
#endif
        return _insert_seq.get();
    }

#if 0
    // insert sequence for all contig and tier1 observations:
    void
    add_insert_seq_obs(//const align_id_t id,
        const std::string& seq) {

        //if(all_read_ids.find(id) != all_read_ids.end()) return;
        _insert_seq.add_obs(seq);
    }

    // dump the full insert_seq_manager
    //
    // this is useful to construct partial indel_data clones
    const insert_seq_manager&
    get_insert_seq_manager() const { return _insert_seq; }
#endif

private:
    // indel key is maintained for debugging only:
    const indel_key _ik;

    mutable insert_seq_manager _insert_seq;

public:

    // map_read_ids refers to the read(s) that support the indel
    // through their genomic alignments, given that those alignments
    // meet the snp-caller's mapping thresholds.
    //
    // submap_read_ids are genomic alignments that fall below the
    // mapping thresholds.
    //
    // all_read_ids contains the list of reads which either have a
    // genomic alignment passing the mapping criteria or have a contig
    // alignment.
    //
    // contig_ids are the grouper contig(s) which support the indel
    //
    typedef std::set<align_id_t> evidence_t;
    evidence_t contig_ids;
    evidence_t all_read_ids; // all contig and tier1 read_ids
    evidence_t tier2_map_read_ids;
    evidence_t submap_read_ids;

    // noise_read_ids indicates that an input alignment had a mismatch
    // fraction which was too high for the alignment to qualify as
    // support for the indel.
    //
    evidence_t noise_read_ids;

    // candidates can be provided from external sources as well:
    bool is_external_candidate;

    // if true candidates should be output even if very unlikely:
    bool is_forced_output;

    // this structure represents support for the indel among all reads
    // which cross an indel breakpoint by a sufficient margin after
    // re-alignment:
    //
    typedef std::map<align_id_t,read_path_scores> score_t;
    score_t read_path_lnp;
    evidence_t suboverlap_tier1_read_ids; // the reads which cross an indel breakpoint, but not by enough
    // to be entered into the scores list
    evidence_t suboverlap_tier2_read_ids;

    ranksum mq_ranksum;
    ranksum baseq_ranksum;
    ranksum read_pos_ranksum;
    unsigned n_mapq;
    // sum of mapq for all reads at this position
    int cumm_mapq;

    struct status_t {
        status_t()
            : is_candidate_indel_cached(false)
            , is_candidate_indel(false)
        {}

        bool is_candidate_indel_cached;
        bool is_candidate_indel;
    };

    mutable status_t status;

private:
#if 0
    // a = a U b
    static
    void
    add_evidence(evidence_t& a,
                 const evidence_t& b) {
        a.insert(b.begin(),b.end());
    }
#endif

    // add observation for the non-shared case
    void
    add_observation_core(const indel_observation_data& obs_data,
                         bool& is_repeat_obs) {

        is_external_candidate=obs_data.is_external_candidate;
        is_forced_output=obs_data.is_forced_output;

        if (! is_external_candidate) {

            using namespace INDEL_ALIGN_TYPE;

            if       (obs_data.iat == CONTIG) {
                contig_ids.insert(obs_data.id);
            } else if (obs_data.is_noise) {
                // noise state overrides all except contig type:
                //
                noise_read_ids.insert(obs_data.id);
            } else if (obs_data.iat == CONTIG_READ) {
                if (all_read_ids.find(obs_data.id) != all_read_ids.end()) {
                    is_repeat_obs=true;
                }
                all_read_ids.insert(obs_data.id);
            } else if (obs_data.iat == GENOME_TIER1_READ) {
                if (all_read_ids.find(obs_data.id) != all_read_ids.end()) {
                    is_repeat_obs=true;
                }
                all_read_ids.insert(obs_data.id);
            } else if (obs_data.iat == GENOME_TIER2_READ) {
                tier2_map_read_ids.insert(obs_data.id);
            } else if (obs_data.iat == GENOME_SUBMAP_READ) {
                submap_read_ids.insert(obs_data.id);
            } else {
                assert(0);
            }
        }
    }

};



// Debugging dumps:
std::ostream& operator<<(std::ostream& os, const indel_observation_data& id);
std::ostream& operator<<(std::ostream& os, const read_path_scores& rps);
std::ostream& operator<<(std::ostream& os, const indel_data& id);


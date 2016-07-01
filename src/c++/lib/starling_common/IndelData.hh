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

#include "blt_util/LogValuePair.hh"
#include "blt_util/ranksum.hh"
#include "starling_common/indel_align_type.hh"
#include "starling_common/IndelKey.hh"
#include "starling_common/starling_types.hh"

#include <cassert>

#include <iosfwd>
#include <map>
#include <string>
#include <set>
#include <vector>


/// represents the data associated with a single indel observation:
///
struct IndelObservationData
{
    bool is_noise = false;
    bool is_external_candidate = false; ///< if true, the allele is automatically promoted to candidate status
    bool is_forced_output = false; ///< if true, the allele must be scored in output
    bool is_discovered_in_active_region = false;
    INDEL_ALIGN_TYPE::index_t iat = INDEL_ALIGN_TYPE::GENOME_SUBMAP_READ;
    align_id_t id = 0;

    std::string breakpointInsertionSequence; ///< partial breakpont sequence, do not use to describe regular indels
};



/// Holds the alignment scores created by each read aligning across an
/// indel which express the relative probability of the read aligning
/// to the indel or the reference (or elsewhere in the genome). Used
/// for indel genotyping.
///
struct ReadPathScores
{
    typedef float score_t;

    ReadPathScores(
        const score_t r=0,
        const score_t i=0,
        const uint16_t ns=0,
        const uint16_t rlen=0,
        const bool is_t1=true,
        const bool is_fwd=true,
        const int16_t rp=0)
        : ref(r)
        , indel(i)
        , nsite(ns)
        , read_length(rlen)
        , is_tier1_read(is_t1)
        , is_fwd_strand(is_fwd)
        , read_pos(rp)
    {}

    void
    insertAlt(
        const IndelKey& indelKey,
        const score_t a);

    score_t ref;
    score_t indel;
    uint16_t nsite;

//    score_t alt;

    // store up to 2 highest scoring alternate indels
    typedef std::vector<std::pair<IndelKey,score_t> > alt_indel_t;
    alt_indel_t alt_indel;

    // used to set expected het allele ratio:
    uint16_t read_length;

    // used to filter for/against tier2 data:
    bool is_tier1_read;

    // so we're able to collect scores by strand
    bool is_fwd_strand;

    // this is used in strelka to calculate read pos ranksums for indels
    int16_t read_pos;
};



/// Accumulates evidence of a consensus breakpoint insert sequence.
///
/// Assumes that long sequences for breakpoints will be noisy and that
/// breakpoints are rare, such that there will be only one insert
/// sequence for a given breakpoint type and position.
///
/// given multiple candidates the "consensus" is just a selection for
/// the longeest insertion, and then a selection of the most common
/// observation among the longest candidates.
///
struct BreakpointInsertSequenceManager
{
    /// get final consensus insert sequence:
    const std::string&
    get()
    {
        if (! _is_consensus)
        {
            _finalize();
        }
        return _consensus_seq;
    }

    // return the insert size
    //
    // note that this test will not trigger consensus finalization,
    // so it is substantially different than get().size()
    unsigned
    getSize() const
    {
        if (_is_consensus) return _consensus_seq.size();

        unsigned size(0);
        for (const auto& val : _obs)
        {
            if (val.first.size() <= size) continue;
            size = val.first.size();
        }

        return size;
    }

    // add insert sequence observation:
    void
    addObservation(const std::string& seq)
    {
        if (seq.empty())
        {
            _exception("Attempting to add breakpoint observation with empty insertion sequence");
        }

        if (_is_consensus)
        {
            _exception("Attempting to add breakpoint insert observation after finalizing");
        }

        // if we don't know what the most common insert will be after
        // this many samples are collected, then more cases will be
        // unlikely to help:
        static const unsigned max_obs_count(256);
        if (_obs_count>max_obs_count) return;

        obs_t::iterator i(_obs.find(seq));
        if (i == _obs.end())
        {
            _obs[seq] = 1;
        }
        else
        {
            i->second += 1;
        }
        _obs_count += 1;
    }

private:
    void _exception(const char* msg) const;

    void _finalize();


    typedef std::map<std::string,unsigned> obs_t;
    bool _is_consensus = false;
    std::string _consensus_seq;
    unsigned _obs_count = 0;
    obs_t _obs;
};


/// indel allele data specific to one sample
struct IndelSampleData
{
    void
    addIndelObservation(
        const IndelObservationData& obs_data);

// ------------- data ------------------

    // tier[12]_map_read ids refers to reads meeting either tier1 or
    // tier2 mapping criteria. All other (non-noise) observations are
    // categorized as submapped
    //
    typedef std::set<align_id_t> evidence_t;
    evidence_t tier1_map_read_ids;
    evidence_t tier2_map_read_ids;
    evidence_t submap_read_ids;

    // noise_read_ids indicates that an input alignment had a mismatch
    // fraction which was too high for the alignment to qualify as
    // support for the indels
    evidence_t noise_read_ids;

    // enumerates support for the indel among all reads
    // which cross an indel breakpoint by a sufficient margin after
    // re-alignment:
    typedef std::map<align_id_t,ReadPathScores> score_t;
    score_t read_path_lnp;

    // the reads which cross an indel breakpoint, but not by enough
    // to be entered into the scores list
    evidence_t suboverlap_tier1_read_ids;
    evidence_t suboverlap_tier2_read_ids;
};


/// indel allele data which is shared across all samples
struct IndelData
{
    IndelData(
        const unsigned sampleCount,
        const IndelKey& indelKey)
        : _sampleData(sampleCount),
          _indelKey(indelKey)
    {
        assert(sampleCount>=1);
    }

    unsigned
    getSampleCount() const
    {
        return _sampleData.size();
    }

    IndelSampleData&
    getSampleData(const unsigned sampleId)
    {
        assert(sampleId<_sampleData.size());
        return _sampleData[sampleId];
    }

    const IndelSampleData&
    getSampleData(const unsigned sampleId) const
    {
        assert(sampleId<_sampleData.size());
        return _sampleData[sampleId];
    }

    void
    addIndelObservation(
        const unsigned sampleId,
        const IndelObservationData& obs_data);

    const std::string&
    getBreakpointInsertSeq() const
    {
        return _breakpointInsertSeq.get();
    }

    /// this test is different than asking for insert-seq in that it does not
    /// trigger insert sequence consensus generation:
    unsigned
    getBreakpointInsertSize() const
    {
        return _breakpointInsertSeq.getSize();
    }

private:
    friend std::ostream& operator<<(std::ostream& os, const IndelData& indelData);
public:
// ------------- data ------------------

    // candidates can be provided from external sources as well:
    bool is_external_candidate = false;

    // if true candidates should be output even if very unlikely:
    bool is_forced_output = false;

    bool isDiscoveredInActiveRegion = false;

    struct ErrorRates
    {
        bool isInit = false;
        LogValuePair indelToRefErrorProb;
        LogValuePair refToIndelErrorProb;

        LogValuePair scaledIndelToRefErrorProb;
        LogValuePair scaledRefToIndelErrorProb;
    };

    struct status_t
    {
        bool is_candidate_indel_cached = false;
        bool is_candidate_indel = false;
    };

    mutable status_t status;
    mutable ErrorRates errorRates;

private:
    std::vector<IndelSampleData> _sampleData;
    mutable BreakpointInsertSequenceManager _breakpointInsertSeq;

    // indel key stored for debugging only
    IndelKey _indelKey;
};



// Debugging dumps:
std::ostream& operator<<(std::ostream& os, const IndelObservationData& id);
std::ostream& operator<<(std::ostream& os, const ReadPathScores& rps);
std::ostream& operator<<(std::ostream& os, const IndelSampleData& indelSampleData);
std::ostream& operator<<(std::ostream& os, const IndelData& indelData);

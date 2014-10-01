// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Starka
// Copyright (c) 2009-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

///
/// \author Chris Saunders
///

#include "blt_util/blt_exception.hh"
#include "starling_common/indel_data.hh"

#include <iostream>
#include <sstream>

//#define DEBUG_ID


void
read_path_scores::insert_alt(const indel_key& ik,
                             const score_t a)
{
    const unsigned ais(static_cast<unsigned>(alt_indel.size()));
    if (ais < 2)
    {
        alt_indel.push_back(std::make_pair(ik,a));
    }
    else
    {
        unsigned min_index(ais);
        score_t min(a);
        for (unsigned i(0); i<ais; ++i)
        {
            if (alt_indel[i].second < min)
            {
                min = alt_indel[i].second;
                min_index = i;
            }
        }
        if (min_index<ais)
        {
            alt_indel[min_index] = std::make_pair(ik,a);
        }
    }
//    log_os << ik << "\n";
}


void
indel_data::add_observation(const indel_observation_data& obs_data,
                            const bool is_shared,
                            bool& is_repeat_obs)
{
#ifdef DEBUG_ID
    log_os << "KATTER: adding obs for indel: " << _ik;
    log_os << "KATTER: is_shared: " << is_shared << " is_repeat: " << is_repeat_obs << "\n";
    log_os << "KATTER: is_external: " << obs_data.is_external_candidate << " align_id: " << obs_data.id << "\n\n";
#endif

    if (! is_shared)
    {
        add_observation_core(obs_data,is_repeat_obs);
    }

    if (! obs_data.insert_seq.empty())
    {
        if (! (is_shared && is_repeat_obs))
        {
            _insert_seq.add_obs(obs_data.insert_seq);
        }
    }
}

// add observation for the non-shared case
void
indel_data::add_observation_core(const indel_observation_data& obs_data,
                                 bool& is_repeat_obs)
{
#ifdef DEBUG_ID
    log_os << "KATTER: adding obs for indel: " << _ik;
    log_os << "KATTER: is_shared: " << is_shared << " is_repeat: " << is_repeat_obs << "\n";
    log_os << "KATTER: is_external: " << obs_data.is_external_candidate << " align_id: " << obs_data.id << "\n\n";
#endif

    // never reset the flags to false if they are true already
    if (! is_external_candidate) is_external_candidate=obs_data.is_external_candidate;
    if (! is_forced_output) is_forced_output=obs_data.is_forced_output;

    if (!is_external_candidate && !is_forced_output)
    {
        using namespace INDEL_ALIGN_TYPE;

        if (obs_data.is_noise)
        {
            // noise state overrides all except contig type:
            //
            noise_read_ids.insert(obs_data.id);
        }
        else if (obs_data.iat == GENOME_TIER1_READ)
        {
            if (all_read_ids.find(obs_data.id) != all_read_ids.end())
            {
                is_repeat_obs=true;
            }
            all_read_ids.insert(obs_data.id);
        }
        else if (obs_data.iat == GENOME_TIER2_READ)
        {
            tier2_map_read_ids.insert(obs_data.id);
        }
        else if (obs_data.iat == GENOME_SUBMAP_READ)
        {
            submap_read_ids.insert(obs_data.id);
        }
        else
        {
            assert(false && "Unknown indel alignment type");
        }
    }
}



std::ostream&
operator<<(std::ostream& os,
           const indel_observation_data& obs)
{

    os << "is_noise: " << obs.is_noise << "\n";
    os << "is_external: " << obs.is_external_candidate << "\n";
    os << "is_forced_output: " << obs.is_forced_output << "\n";
    os << "type: " << INDEL_ALIGN_TYPE::label(obs.iat) << "\n";
    os << "align_id: " << obs.id << "\n";
    os << "insert_seq: " << obs.insert_seq << "\n";
    return os;
}

static
void
report_indel_evidence_set(const indel_data::evidence_t& e,
                          const char* label,
                          std::ostream& os)
{
    typedef indel_data::evidence_t::const_iterator viter;
    viter i(e.begin()),i_end(e.end());
    for (unsigned n(0); i!=i_end; ++i)
    {
        os << label << " no: " << ++n << " id: " << *i << "\n";
    }
}





std::ostream&
operator<<(std::ostream& os,
           const read_path_scores& rps)
{

    os << "ref: " << rps.ref
       << " indel: " << rps.indel
       << " nsite: " << rps.nsite;

#if 0
    if (rps.is_alt)
    {
        os << " alt: " << rps.alt;
    }
#else
    typedef read_path_scores::alt_indel_t::const_iterator aiter;
    aiter i(rps.alt_indel.begin()), i_end(rps.alt_indel.end());
    for (; i!=i_end; ++i)
    {
        const indel_key& ik(i->first);
        os << " alt-" << ik.pos << "-" << INDEL::get_index_label(ik.type) << ik.length << ": " << i->second;
    }
#endif

    os << " ist1?: " << rps.is_tier1_read;

    return os;
}


void
insert_seq_manager::
_exception(const char* msg) const
{
    std::ostringstream oss;
    oss << "Exception in insert_seq_manager: " << msg;
    throw blt_exception(oss.str().c_str());
}


void
insert_seq_manager::
_finalize()
{
    unsigned count(0);
    std::string& candidate(_consensus_seq);

    for (const auto& val : _obs)
    {
        if (val.first.size() < candidate.size()) continue;
        if (val.first.size() == candidate.size())
        {
            if (val.second <= count) continue;
        }
        candidate = val.first;
        count = val.second;
    }
    _consensus_seq = candidate;
    _obs.clear();
    _is_consensus=true;
}



std::ostream&
operator<<(std::ostream& os,
           const indel_data& id)
{
    os << "seq: " << id.get_insert_seq() << "\n";

    report_indel_evidence_set(id.all_read_ids,"all_read",os);
    //    report_indel_evidence_set(id.tier1_map_read_ids,"tier1_map_read",os);
    report_indel_evidence_set(id.tier2_map_read_ids,"tier2_map_read",os);
    report_indel_evidence_set(id.submap_read_ids,"submap_read",os);
    report_indel_evidence_set(id.noise_read_ids,"noise_read",os);

    {
        typedef indel_data::score_t::const_iterator siter;
        siter i(id.read_path_lnp.begin()), i_end(id.read_path_lnp.end());
        for (unsigned n(0); i!=i_end; ++i)
        {
            os << "read_path_lnp no: " << ++n
               << " id: " << i->first
               << " " << i->second
               << "\n";
        }
    }

    report_indel_evidence_set(id.suboverlap_tier1_read_ids,"suboverlap_tier1_read",os);
    report_indel_evidence_set(id.suboverlap_tier2_read_ids,"suboverlap_tier2_read",os);

    os << "is_external_candidate: " << id.is_external_candidate << "\n";

    return os;
}

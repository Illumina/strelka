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

#include "SequenceErrorCountsPosProcessor.hh"
#include "blt_common/ref_context.hh"
#include "common/OutStream.hh"
#include "starling_common/starling_indel_call_pprob_digt.hh"
#include "starling_common/starling_diploid_indel.hh"
#include "starling_common/starling_indel_error_prob.hh"

#include <algorithm>
#include <vector>



SequenceErrorCountsPosProcessor::
SequenceErrorCountsPosProcessor(
    const SequenceErrorCountsOptions& opt,
    const SequenceErrorCountsDerivOptions& dopt,
    const reference_contig_segment& ref,
    const SequenceErrorCountsStreams& streams)
    : base_t(opt,dopt,ref,streams,1),
      _opt(opt),
      _streams(streams)
{
    // check that we have write permission on the output file early:
    {
        OutStream outs(opt.countsFilename);
        if (opt.is_write_observations())
        {
            OutStream outs2(opt.observationsBedFilename);
        }
    }

    // setup indel syncronizers:
    {
        sample_info& normal_sif(sample(0));

        if (dopt.is_max_depth())
        {
            if (opt.max_candidate_indel_depth_factor > 0.)
            {
                _max_candidate_normal_sample_depth = (opt.max_candidate_indel_depth_factor * dopt.max_depth);
            }
        }

        if (opt.max_candidate_indel_depth > 0.)
        {
            if (_max_candidate_normal_sample_depth > 0.)
            {
                _max_candidate_normal_sample_depth = std::min(_max_candidate_normal_sample_depth,static_cast<double>(opt.max_candidate_indel_depth));
            }
            else
            {
                _max_candidate_normal_sample_depth = opt.max_candidate_indel_depth;
            }
        }

        indel_sync_data isdata;
        isdata.register_sample(normal_sif.indel_buff,normal_sif.estdepth_buff,normal_sif.estdepth_buff_tier2,
                               normal_sif.sample_opt, _max_candidate_normal_sample_depth, 0);
        normal_sif.indel_sync_ptr.reset(new indel_synchronizer(opt, ref, dopt.countCache, isdata, 0));
    }
}

SequenceErrorCountsPosProcessor::
~SequenceErrorCountsPosProcessor()
{
    _counts.save(_opt.countsFilename.c_str());
}



void
SequenceErrorCountsPosProcessor::
reset()
{
    base_t::reset();
}



/// generalization of overlapping indels:
///
struct OrthogonalHaplotypeCandidateGroup
{
    const indel_key&
    key(
        const unsigned index) const
    {
        return iter(index)->first;
    }

    const indel_data&
    data(
        const unsigned index) const
    {
        assert(index < size());
        return get_indel_data(iter(index));
    }

    const indel_buffer::const_iterator&
    iter(
        const unsigned index) const
    {
        assert(index < size());
        return haps[index];
    }

    unsigned
    size() const
    {
        return haps.size();
    }

    bool
    empty() const
    {
        return haps.empty();
    }

    void
    clear()
    {
        haps.clear();
    }

    void
    addHaplotype(
        const indel_buffer::const_iterator hapIter)
    {
        haps.push_back(hapIter);
    }

    std::vector<indel_buffer::const_iterator> haps;
};




#include "blt_util/math_util.hh"
#include "blt_util/prob_util.hh"

/// order overlapping indels by total read support
///
/// \param[out] referenceRank rank of the reference among all haplotypes, with 0 being the best (highest scoring)
///
static
void
getOrthogonalHaplotypeSupportCounts(
    const OrthogonalHaplotypeCandidateGroup& hg,
    std::vector<unsigned>& support)
{
    const unsigned nonrefHapCount(hg.size());
    assert(nonrefHapCount!=0);

    std::set<unsigned> readIds;

    // intersection of read ids
    {
        std::map<unsigned,unsigned> countReadIds;
        for (unsigned nonrefHapIndex(0); nonrefHapIndex<nonrefHapCount; nonrefHapIndex++)
        {
            const indel_data& id(hg.data(nonrefHapIndex));
            for (const auto& score : id.read_path_lnp)
            {
                const auto iter(countReadIds.find(score.first));
                if (iter==countReadIds.end())
                {
                    countReadIds.insert(std::make_pair(score.first,1));
                }
                else
                {
                    iter->second += 1;
                }
            }
        }

        for (const auto& value : countReadIds)
        {
            if (value.second >= nonrefHapCount)
            {
                readIds.insert(value.first);
            }
        }
    }


    // count of all haplotypes including reference
    const unsigned allHapCount(nonrefHapCount+1);
    const unsigned refHapIndex(nonrefHapCount);

    support.clear();
    support.resize(allHapCount,0);

    std::vector<double> lhood(allHapCount);
    for (const auto readId : readIds)
    {
        static const double log0(-std::numeric_limits<double>::infinity());
        std::fill(lhood.begin(),lhood.end(),log0);
        for (unsigned nonrefHapIndex(0); nonrefHapIndex<nonrefHapCount; nonrefHapIndex++)
        {
            const indel_data& id(hg.data(nonrefHapIndex));

            const auto iditer(id.read_path_lnp.find(readId));
            if (iditer==id.read_path_lnp.end())
            {
                continue;
            }
            const read_path_scores& path_lnp(iditer->second);

            lhood[refHapIndex] = std::max(lhood[refHapIndex],static_cast<double>(path_lnp.ref));
            lhood[nonrefHapIndex] = path_lnp.indel;
        }
        unsigned maxIndex(0);
        normalize_ln_distro(lhood.begin(),lhood.end(),maxIndex);

        for (unsigned allHapIndex(0); allHapIndex<allHapCount; allHapIndex++)
        {
            if (lhood[allHapIndex]>0.999)
            {
                support[allHapIndex]++;
            }
        }
    }
}


static
INDEL_SIGNAL_TYPE::index_t
getIndelType(
    const starling_indel_report_info& iri)
{
    int rudiff(static_cast<int>(iri.ref_repeat_count)-static_cast<int>(iri.indel_repeat_count));
    rudiff = std::min(3,std::max(-3,rudiff));

    using namespace INDEL_SIGNAL_TYPE;
    switch (rudiff)
    {
    case  1:
        return DELETE_1;
    case  2:
        return DELETE_2;
    case  3:
        return DELETE_GE3;
    case -1:
        return INSERT_1;
    case -2:
        return INSERT_2;
    case -3:
        return INSERT_GE3;
    default:
        assert(false);
        return SIZE;
    }
}



static
void
mergeIndelObservations(
    const IndelErrorContext& context,
    const IndelErrorContextObservation& obs,
    std::map<IndelErrorContext,IndelErrorContextObservation>& indelObservations)
{
    auto iter(indelObservations.find(context));

    if (iter == indelObservations.end())
    {
        indelObservations.insert(std::make_pair(context,obs));
    }
    else
    {
        // all signal counts should be summed:
        for (unsigned signalIndex(0); signalIndex<INDEL_SIGNAL_TYPE::SIZE; ++signalIndex)
        {
            iter->second.signalCounts[signalIndex] += obs.signalCounts[signalIndex];
        }

        // refCount used for the group is the lowest submitted:
        iter->second.refCount = std::min(iter->second.refCount, obs.refCount);
    }
}



void
SequenceErrorCountsPosProcessor::
process_pos_error_counts(
    const pos_t pos,
    const unsigned sample_no)
{
    static const unsigned maxHpolLength(20);

    // note multi-sample status -- can still be called only for one sample
    // and only for sample 0. working on generalization:
    //
    if (sample_no!=0) return;

    const char refBase(_ref.get_base(pos));

    if (refBase=='N') return;

    IndelErrorCounts& indelCounts(_counts.getIndelCounts());

    // Current multiploid indel model can handle a het or hom indel
    // allele vs. reference, or two intersecting non-reference indel
    // alleles. (note that indel intersection is evaluated only in
    // terms of breakpoints -- so, for instance, a small het deletion
    // could occur within a large het deletion and the two would be
    // treated as non-interacting -- this is just an artifact of how
    // the methods are coded,)
    //

    typedef indel_buffer::const_iterator ciiter;

    sample_info& sif(sample(sample_no));


    // candidate indels are turned off above a certain depth relative to chromosome
    // mean, so we make sure to turn off counts as well to prevent a reference bias:
    {
        if (_max_candidate_normal_sample_depth > 0.)
        {
            const unsigned estdepth(sif.estdepth_buff.val(pos-1));
            const unsigned estdepth2(sif.estdepth_buff_tier2.val(pos-1));
            if ((estdepth+estdepth2) > _max_candidate_normal_sample_depth)
            {
                // record the number of contexts which have been skipped due to high depth:
                IndelErrorContext context;
                context.repeatCount = 1;
                indelCounts.addDepthSkip(context);

                const unsigned leftHpolSize(get_left_shifted_hpol_size(pos,_ref));
                if (leftHpolSize>1)
                {
                    context.repeatCount = std::min(maxHpolLength,leftHpolSize);
                    indelCounts.addDepthSkip(context);
                }
                return;
            }
        }
    }

    // define groups of overlapping indels:
    //
    // indels which form "conflict cliques" go into a single indel
    // group for scoring together
    //
    // this design should look overblown within this function, it
    // is intended to generalize across multiple positions, in which
    // case we will have to deal with non-clique conflict patterns.
    //
    ciiter it(sif.indel_buff.pos_iter(pos));
    const ciiter it_end(sif.indel_buff.pos_iter(pos+1));

    std::vector<OrthogonalHaplotypeCandidateGroup> _groups;

    for (; it!=it_end; ++it)
    {
        const indel_key& ik(it->first);
        const indel_data& id(get_indel_data(it));

        if (ik.is_breakpoint()) continue;

        const bool isForcedOutput(id.is_forced_output);

        if (! isForcedOutput)
        {
            if (! sif.indel_sync().is_candidate_indel(ik,id)) continue;
        }

        if (! (ik.type == INDEL::DELETE || ik.type == INDEL::INSERT)) continue;

        // all indels at the same position are conflicting:
        if (_groups.empty()) _groups.resize(1);
        OrthogonalHaplotypeCandidateGroup& hg(_groups[0]);
        hg.addHaplotype(it);
    }

    // buffer observations until we get through all overlaping indels at this site:
    std::map<IndelErrorContext,IndelErrorContextObservation> indelObservations;


    if (! _groups.empty())
    {
        const OrthogonalHaplotypeCandidateGroup& hg(_groups[0]);
        const unsigned nonrefHapCount(hg.size());
        std::vector<unsigned> support;
        getOrthogonalHaplotypeSupportCounts(hg,support);

        for (unsigned nonrefHapIndex(0); nonrefHapIndex<nonrefHapCount; ++nonrefHapIndex)
        {
            const indel_key& ik(hg.key(nonrefHapIndex));
            const indel_data& id(hg.data(nonrefHapIndex));

            starling_indel_report_info iri;
            get_starling_indel_report_info(ik,id,_ref,iri);

            IndelErrorContext context;
            if ((iri.repeat_unit_length==1) && (iri.ref_repeat_count>1))
            {
                // guard against the occasional non-normalized indel:
                const unsigned leftHpolSize(get_left_shifted_hpol_size(pos,_ref));
                if (leftHpolSize == iri.ref_repeat_count)
                {
                    context.repeatCount = std::min(maxHpolLength, iri.ref_repeat_count);
                }
            }

            IndelErrorContextObservation obs;

            const INDEL_SIGNAL_TYPE::index_t sigIndex(getIndelType(iri));
            obs.signalCounts[sigIndex] = support[nonrefHapIndex];
            obs.refCount = support[nonrefHapCount];

            // an indel candidate can have 0 q30 indel reads when it is only supported by
            // noise reads (i.e. indel occurs outside of a read's valid alignment range,
            // see lib/starling_common/starling_read_util.cpp::get_valid_alignment_range)
            // in this case, we're not going to report the incidence as noise, since it's
            // not a read we would consider in variant calling
            if (support[nonrefHapIndex] > 0 && _opt.is_write_observations())
            {
                std::ostream& obs_os(*_streams.observation_bed_osptr());
                obs_os << _opt.bam_seq_name << "\t";
                obs_os << ik.pos << "\t" << ik.pos + ik.length << "\t" << INDEL::get_index_label(ik.type) << "\t";
                obs_os << iri.repeat_unit << "\t" << iri.ref_repeat_count << "\t";
                obs_os << ik.length << "\t" << nonrefHapIndex << "/" << nonrefHapCount << "\t";
                obs_os << support[nonrefHapIndex] << "\t" << support[nonrefHapCount] << "\t";
                obs_os << std::accumulate(support.begin(), support.end(), 0) << std::endl;
            }

            mergeIndelObservations(context,obs,indelObservations);
        }
    }

    // background depth is always one minus position to be consistent with indel report:
    const pos_t depth_pos(pos-1);
    const snp_pos_info& spi(sif.bc_buff.get_pos(depth_pos));
    const unsigned depth(spi.calls.size());

    for (const auto& value : indelObservations)
    {
        indelCounts.addError(value.first,value.second,depth);
    }

    // add all the backgrounds that haven't been covered already
    {
        unsigned leftHpolSize(get_left_shifted_hpol_size(pos,_ref));

        IndelErrorContext context;
        IndelBackgroundObservation obs;
        obs.depth = depth;

        // always add hpol=1:
        if (! indelObservations.count(context))
        {
            indelCounts.addBackground(context,obs);
        }

        // also check for a more specific reference context
        if (leftHpolSize>1)
        {
            context.repeatCount = std::min(maxHpolLength,leftHpolSize);
            if (! indelObservations.count(context))
            {
                indelCounts.addBackground(context,obs);
            }
        }
    }
}

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

    static const unsigned sampleId(0);

    // setup indel syncronizer:
    {
        sample_info& normal_sif(sample(sampleId));

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

        sample_id_t syncSampleId;
        syncSampleId = indel_sync().register_sample(normal_sif.estdepth_buff, normal_sif.estdepth_buff_tier2,
                                                    normal_sif.sample_opt, _max_candidate_normal_sample_depth);

        assert(syncSampleId == sampleId);

        indel_sync().finalizeSamples();
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



void
SequenceErrorCountsPosProcessor::
insertExcludedRegion(
    const known_pos_range2& excludedRange)
{
    _stageman.validate_new_pos_value(excludedRange.begin_pos(),STAGE::READ_BUFFER);
    _excludedRegions.addRegion(excludedRange);
    _is_skip_process_pos=false;
}

void
SequenceErrorCountsPosProcessor::
addKnownVariant(
    const vcf_record& knownVariant)
{
    _stageman.validate_new_pos_value(knownVariant.pos, STAGE::READ_BUFFER);
    _knownVariants.addVcfRecord(knownVariant);
    _is_skip_process_pos=false;
}



static
bool
isKnownVariantMatch(
    const RecordTracker::indel_value_t& knownVariants,
    const indel_key& ik,
    const IndelData& id,
    RecordTracker::indel_value_t& overlap)
{
    if (knownVariants.empty()) return false;

    // the variant matches a known variant
    // if it overlaps a known record and has
    // the same type and length.  If the
    // known variant is an insertion, then
    // the inserted sequence must match
    for (const auto& kv : knownVariants)
    {
        if (kv.altMatch(ik, id)) overlap.insert(kv);
    }
    return (! overlap.empty());
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

    const IndelData&
    data(
        const unsigned index) const
    {
        assert(index < size());
        return get_indel_data(iter(index));
    }

    const indel_synchronizer::const_iterator&
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
        const indel_synchronizer::const_iterator hapIter)
    {
        haps.push_back(hapIter);
    }

    std::vector<indel_synchronizer::const_iterator> haps;
};




#include "blt_util/math_util.hh"
#include "blt_util/prob_util.hh"

/// order overlapping indels by total read support in one sample
///
/// \param[out] referenceRank rank of the reference among all haplotypes, with 0 being the best (highest scoring)
///
static
void
getOrthogonalHaplotypeSupportCounts(
    const OrthogonalHaplotypeCandidateGroup& hg,
    const unsigned sampleId,
    std::vector<unsigned>& support,
    const bool is_tier1_only = true)
{
    const unsigned nonrefHapCount(hg.size());
    assert(nonrefHapCount!=0);

    std::set<unsigned> readIds;

    // intersection of read ids which have a likelihood evaluated over all candidate haplotypes:
    {
        std::map<unsigned,unsigned> countReadIds;
        for (unsigned nonrefHapIndex(0); nonrefHapIndex<nonrefHapCount; nonrefHapIndex++)
        {
            const IndelSampleData& isd(hg.data(nonrefHapIndex).getSampleData(sampleId));
            for (const auto& score : isd.read_path_lnp)
            {
                if (is_tier1_only && (! score.second.is_tier1_read)) continue;

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
            const IndelSampleData& isd(hg.data(nonrefHapIndex).getSampleData(sampleId));

            const auto iditer(isd.read_path_lnp.find(readId));

            if (iditer==isd.read_path_lnp.end())
            {
                // note we've set readIds to intersection now so this should never happen:
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

        // known variant status is the maximum status within the set (e.g. if UNKNOWN
        // and VARIANT, site is variant)

        iter->second.variantStatus = std::max(iter->second.variantStatus, obs.variantStatus);
    }
}



void
SequenceErrorCountsPosProcessor::
process_pos_error_counts(
    const pos_t pos,
    const unsigned sampleId)
{
    static const unsigned maxHpolLength(20);

    // note multi-sample status -- can still be called only for one sample
    // and only for sample 0. working on generalization:
    //
    if (sampleId!=0) return;

    const char refBase(_ref.get_base(pos));

    if (refBase=='N') return;

    BaseErrorCounts& baseCounts(_counts.getBaseCounts());
    IndelErrorCounts& indelCounts(_counts.getIndelCounts());


    const sample_info& sif(sample(sampleId));

    // right now there's only one baseContext, so just set it const here and leave it
    const BaseErrorContext baseContext;

    bool isSkipSNV(false);
    bool isSkipIndel(false);

    if (_excludedRegions.isIntersectRegion(pos))
    {
        baseCounts.addExcludedRegionSkip(baseContext);
        isSkipSNV=true;

        IndelErrorContext indelContext;
        indelCounts.addExcludedRegionSkip(indelContext);
        const unsigned leftHpolSize(get_left_shifted_hpol_size(pos,_ref));
        if (leftHpolSize>1)
        {
            indelContext.repeatCount = std::min(maxHpolLength,leftHpolSize);
            indelCounts.addExcludedRegionSkip(indelContext);
        }
        isSkipIndel=true;
    }
    else if (_max_candidate_normal_sample_depth > 0.)
    {
        //
        // candidate variants are turned off above a certain depth relative to chromosome mean
        // for basecalls we check the current position, for indels we check the previous position.
        //
        {
            // handle SNVs
            const unsigned estdepth(sif.estdepth_buff.val(pos));
            const unsigned estdepth2(sif.estdepth_buff_tier2.val(pos));
            if ((estdepth+estdepth2) > _max_candidate_normal_sample_depth)
            {
                baseCounts.addDepthSkip(baseContext);
                isSkipSNV=true;
            }
        }

        {
            const unsigned estdepth(sif.estdepth_buff.val(pos-1));
            const unsigned estdepth2(sif.estdepth_buff_tier2.val(pos-1));
            if ((estdepth+estdepth2) > _max_candidate_normal_sample_depth)
            {
                // handle indels:
                IndelErrorContext indelContext;
                indelCounts.addDepthSkip(indelContext);

                const unsigned leftHpolSize(get_left_shifted_hpol_size(pos,_ref));
                if (leftHpolSize>1)
                {
                    indelContext.repeatCount = std::min(maxHpolLength,leftHpolSize);
                    indelCounts.addDepthSkip(indelContext);
                }
                isSkipIndel=true;
            }
        }
    }


    // handle base error signal
    if (! isSkipSNV)
    {
        // be relatively intolerant of anything interesting happening in the local sequence neighborhood:
        static const double snvMMDRMaxFrac(0.05);

        const snp_pos_info& sinfo(sif.bc_buff.get_pos(pos));
        unsigned fcount(0);
        for (const base_call& bc : sinfo.calls)
        {
            if (bc.is_call_filter) fcount++;
        }
        fcount += sinfo.tier2_calls.size();

        const double filtFrac(fcount/static_cast<double>(sinfo.calls.size()+sinfo.tier2_calls.size()));

        if (filtFrac >= snvMMDRMaxFrac)
        {
            baseCounts.addNoiseSkip(baseContext);
        }
        else
        {
            bool isEmpty(true);
            const uint8_t ref_id(base_to_id(refBase));
            BaseErrorContextInputObservation obs;
            for (const base_call& bc : sinfo.calls)
            {
                if (bc.is_call_filter) continue;

                const uint16_t qual(bc.get_qscore());
                if (qual<_opt.min_qscore) continue;

                static const uint16_t min_count_qscore(25);
                if (qual<min_count_qscore) continue;
                isEmpty=false;
                const bool isRef(bc.base_id == ref_id);
                if (isRef)
                {
                    obs.addRefCount(bc.is_fwd_strand, bc.get_qscore());
                }
                else
                {
                    obs.addAltCount(bc.is_fwd_strand, bc.get_qscore());
                }
            }

            if (isEmpty)
            {
                baseCounts.addEmptySkip(baseContext);
            }
            else
            {
                baseCounts.addSiteObservation(baseContext,obs);
            }
        }
    }


    if (isSkipIndel) return;


    // define groups of overlapping indels:
    //
    // indels which form "conflict cliques" go into a single indel
    // group for scoring together
    //
    // this design should look overblown within this function, it
    // is intended to generalize across multiple positions, in which
    // case we will have to deal with non-clique conflict patterns.
    //
    std::vector<OrthogonalHaplotypeCandidateGroup> _groups;

    auto it(indel_sync().pos_iter(pos));
    const auto it_end(indel_sync().pos_iter(pos+1));

    for (; it!=it_end; ++it)
    {
        const indel_key& ik(it->first);
        const IndelData& id(get_indel_data(it));

        if (ik.is_breakpoint()) continue;

        const bool isForcedOutput(id.is_forced_output);

        if (! isForcedOutput)
        {
            if (! indel_sync().is_candidate_indel(ik, id)) continue;
        }

        if (! (ik.type == INDEL::DELETE || ik.type == INDEL::INSERT)) continue;

        // all indels at the same position are conflicting:
        if (_groups.empty()) _groups.resize(1);
        OrthogonalHaplotypeCandidateGroup& hg(_groups[0]);
        hg.addHaplotype(it);
    }

    // buffer observations until we get through all overlaping indels at this site:
    std::map<IndelErrorContext,IndelErrorContextObservation> indelObservations;


    // check for any known variants overlapping this position
    RecordTracker::indel_value_t knownVariantRecords;
    _knownVariants.intersectingRecord(pos, knownVariantRecords);

    if (! _groups.empty())
    {
        const OrthogonalHaplotypeCandidateGroup& hg(_groups[0]);
        const unsigned nonrefHapCount(hg.size());
        std::vector<unsigned> support;
        getOrthogonalHaplotypeSupportCounts(hg, sampleId, support);

        for (unsigned nonrefHapIndex(0); nonrefHapIndex<nonrefHapCount; ++nonrefHapIndex)
        {
            const indel_key& ik(hg.key(nonrefHapIndex));
            const IndelData& id(hg.data(nonrefHapIndex));

            starling_indel_report_info iri;
            get_starling_indel_report_info(ik, id,_ref,iri);

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

            RecordTracker::indel_value_t overlappingRecords;
            isKnownVariantMatch(knownVariantRecords, ik, id, overlappingRecords);

            IndelErrorContextObservation obs;

            const INDEL_SIGNAL_TYPE::index_t sigIndex(getIndelType(iri));
            obs.signalCounts[sigIndex] = support[nonrefHapIndex];
            obs.refCount = support[nonrefHapCount];
            obs.assignKnownStatus(overlappingRecords);

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
                obs_os << GENOTYPE_STATUS::label(obs.variantStatus) << "\t";
                obs_os << context.repeatCount << "\t" << INDEL_SIGNAL_TYPE::label(sigIndex) << "\t";
                obs_os << ik.length << "\t" << nonrefHapIndex + 1 << "/" << nonrefHapCount << "\t";
                obs_os << obs.signalCounts[sigIndex] << "\t" << obs.refCount << "\t";
                obs_os << std::accumulate(support.begin(), support.end(), 0) << std::endl;
            }

#if 0
            if (obs.variantStatus == GENOTYPE_STATUS::VARIANT &&
                obs.signalCounts[sigIndex] == 0)
            {
                std::ostream& obs_os(std::cout);
                obs_os << _opt.bam_seq_name << "\t";
                obs_os << ik.pos << "\t" << ik.pos + ik.length << "\t" << INDEL::get_index_label(ik.type) << "\t";
                obs_os << iri.repeat_unit << "\t" << iri.ref_repeat_count << "\t";
                obs_os << GENOTYPE_STATUS::label(obs.variantStatus) << "\t";
                obs_os << ik.length << "\t" << nonrefHapIndex + 1 << "/" << nonrefHapCount << "\t";
                obs_os << support[nonrefHapCount] << "\t";
                obs_os << std::accumulate(support.begin(), support.end(), 0) << std::endl;

                for (const auto& rec : overlappingRecords)
                {
                    std::cout << rec << std::endl;
                }
            }
#endif

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
        obs.assignKnownStatus(knownVariantRecords);
        // the assumption is that a background position should have
        // the variant status of any overlapping known variants,
        // regardless of whether the genotypes match

        // std::cout << "POS: " << pos << std::endl;
        // std::cout << obs << std::endl;

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

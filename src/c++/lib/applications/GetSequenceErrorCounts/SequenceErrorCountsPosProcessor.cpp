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
#include "starling_common/OrthogonalVariantAlleleCandidateGroupUtil.hh"
#include "starling_common/starling_indel_call_pprob_digt.hh"

#include <algorithm>
#include <vector>



SequenceErrorCountsPosProcessor::
SequenceErrorCountsPosProcessor(
    const SequenceErrorCountsOptions& opt,
    const SequenceErrorCountsDerivOptions& dopt,
    const reference_contig_segment& ref,
    const SequenceErrorCountsStreams& streams)
    : base_t(opt,dopt,ref,streams, opt.alignFileOpt.alignmentFilename.size()),
      _opt(opt),
      _streams(streams)
{
    // not generalized to multi-sample yet:
    assert(getSampleCount()==1);
    static const unsigned sampleId(0);

    // check that we have write permission on the output file early:
    {
        OutStream outs(opt.countsFilename);
        if (opt.is_write_observations())
        {
            OutStream outs2(opt.observationsBedFilename);
        }
    }

    // setup indel buffer
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

        getIndelBuffer().setMaxCandidateDepth(_max_candidate_normal_sample_depth);

        const unsigned syncSampleId = getIndelBuffer().registerSample(normal_sif.estdepth_buff, normal_sif.estdepth_buff_tier2, true);

        assert(syncSampleId == sampleId);

        getIndelBuffer().finalizeSamples();
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
    const IndelKey& indelKey,
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
        if (kv.altMatch(indelKey)) overlap.insert(kv);
    }
    return (! overlap.empty());
}



static
void
getOrthogonalHaplotypeSupportCounts(
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    const unsigned sampleId,
    std::vector<unsigned>& support,
    const bool isTier1Only = true)
{
    const unsigned nonrefAlleleCount(alleleGroup.size());
    assert(nonrefAlleleCount!=0);

    // intersection of read ids which have a likelihood evaluated over all candidate haplotypes:
    std::set<unsigned> readIds;
    getAlleleGroupIntersectionReadIds(sampleId, alleleGroup, readIds, isTier1Only);

    // count of all haplotypes including reference
    const unsigned fullAlleleCount(nonrefAlleleCount+1);

    support.clear();
    support.resize(fullAlleleCount,0);

    for (const auto readId : readIds)
    {
        std::vector<double> lhood(fullAlleleCount);
        getAlleleLikelihoodsFromRead(sampleId, alleleGroup, readId, lhood);
        for (unsigned fullAlleleIndex(0); fullAlleleIndex<fullAlleleCount; fullAlleleIndex++)
        {
            if (lhood[fullAlleleIndex]>0.999)
            {
                support[fullAlleleIndex]++;
            }
        }
    }
}



static
INDEL_SIGNAL_TYPE::index_t
getIndelType(
    const starling_indel_report_info& indelReportInfo)
{
    int rudiff(static_cast<int>(indelReportInfo.ref_repeat_count)-static_cast<int>(indelReportInfo.indel_repeat_count));
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


    // define groups of overlapping alleles to rank and then genotype.
    //
    // overlapping alleles can be thought to form "conflict graphs", where an edge exists between two alleles
    // which cannot exist together on the same haplotype (called orthogonal alleles below). Without phasing
    // information, we can only (accurately) genotype among sets of alleles forming a clique in the graph.
    //
    // Given above constraint, we first identify all candidates alleles with a start position at the current
    // allele genotyper position (these form a clique by definition), and then greedily add the top-ranking
    // overlapping alleles with different start positions if they preserve the orthogonal clique relationship
    // of the set.
    //
    // Once we have the largest possible allele set, we skip this locus if the total number of alleles is too
    // large. Otherwise, read support for each allele in the overlapping set is enumberated. For all alleles that
    // start at this position, we record the supporing counts in the indel error analysis counting structure.
    //

    auto it(getIndelBuffer().positionIterator(pos));
    const auto it_end(getIndelBuffer().positionIterator(pos + 1));

    OrthogonalVariantAlleleCandidateGroup orthogonalVariantAlleles;
    for (; it!=it_end; ++it)
    {
        const IndelKey& indelKey(it->first);
        const IndelData& indelData(getIndelData(it));

        if (indelKey.is_breakpoint()) continue;

        const bool isForcedOutput(indelData.isForcedOutput);
        if (not isForcedOutput)
        {
            const IndelSampleData& indelSampleData(indelData.getSampleData(sampleId));
            const bool isZeroCoverage(indelSampleData.read_path_lnp.empty());

            if (isZeroCoverage) continue;
            if (not getIndelBuffer().isCandidateIndel(indelKey, indelData)) continue;
        }

        if (! (indelKey.isPrimitiveDeletionAllele() || indelKey.isPrimitiveInsertionAllele())) continue;

        // all alleles at the same position are automatically conflicting/orthogonal:
        orthogonalVariantAlleles.addVariantAllele(it);
    }


    // this takes the place of the ploidy argument we use during real variant calling -- the
    // question is: should we set this to 2 and filter out the items we normally filter out during
    // variant calling? Or do we want ot try to capture something closer to the true underlying noise
    // distribution?
    //
    // in this case we choose the latter in principal (favor true noise vs what the germline caller seees),
    // but still set an upper-limit on the total number of overlapping variants, recognizing that we can't
    // handle the noise accurately as variant density goes up
    //
    const unsigned maxOverlap(4);


    if (orthogonalVariantAlleles.size() > maxOverlap) return;

    // check for any known variants overlapping this position
    RecordTracker::indel_value_t knownVariantRecords;
    _knownVariants.intersectingRecord(pos, knownVariantRecords);

    // buffer observations until we get through all overlapping indels at this position:
    std::map<IndelErrorContext,IndelErrorContextObservation> indelObservations;

    if (not orthogonalVariantAlleles.empty())
    {
        {
            const bool isEveryAltIncluded = \
                                            addAllelesAtOtherPositions(pos, get_largest_total_indel_ref_span_per_read(), sampleId,
                                                                       (maxOverlap + 1), getIndelBuffer(), orthogonalVariantAlleles);

            if (not isEveryAltIncluded) return;

            if (orthogonalVariantAlleles.size() > maxOverlap) return;
        }

        const unsigned nonrefAlleleCount(orthogonalVariantAlleles.size());
        std::vector<unsigned> support;
        getOrthogonalHaplotypeSupportCounts(orthogonalVariantAlleles, sampleId, support);

        for (unsigned nonrefAlleleIndex(0); nonrefAlleleIndex<nonrefAlleleCount; ++nonrefAlleleIndex)
        {
            const IndelKey& indelKey(orthogonalVariantAlleles.key(nonrefAlleleIndex));

            if (indelKey.pos != pos) continue;

            const IndelData& indelData(orthogonalVariantAlleles.data(nonrefAlleleIndex));

            starling_indel_report_info indelReportInfo;
            get_starling_indel_report_info(indelKey, indelData,_ref,indelReportInfo);

            IndelErrorContext context;
            if ((indelReportInfo.repeat_unit_length==1) && (indelReportInfo.ref_repeat_count>1))
            {
                // guard against the occasional non-normalized indel:
                const unsigned leftHpolSize(get_left_shifted_hpol_size(pos,_ref));
                if (leftHpolSize == indelReportInfo.ref_repeat_count)
                {
                    context.repeatCount = std::min(maxHpolLength, indelReportInfo.ref_repeat_count);
                }
            }

            RecordTracker::indel_value_t overlappingRecords;
            isKnownVariantMatch(knownVariantRecords, indelKey, overlappingRecords);

            IndelErrorContextObservation obs;

            const INDEL_SIGNAL_TYPE::index_t sigIndex(getIndelType(indelReportInfo));
            obs.signalCounts[sigIndex] = support[nonrefAlleleIndex];
            obs.refCount = support[nonrefAlleleCount];
            obs.assignKnownStatus(overlappingRecords);

            // an indel candidate can have 0 q30 indel reads when it is only supported by
            // noise reads (i.e. indel occurs outside of a read's valid alignment range,
            // see lib/starling_common/starling_read_util.cpp::get_valid_alignment_range)
            // in this case, we're not going to report the incidence as noise, since it's
            // not a read we would consider in variant calling
            if (support[nonrefAlleleIndex] > 0 && _opt.is_write_observations())
            {
                std::ostream& obs_os(*_streams.observation_bed_osptr());
                obs_os << _opt.bam_seq_name << "\t";
                obs_os << indelKey.pos << "\t" << indelKey.pos + indelKey.deletionLength << "\t" << INDEL::get_index_label(indelKey.type) << "\t";
                obs_os << indelReportInfo.repeat_unit << "\t" << indelReportInfo.ref_repeat_count << "\t";
                obs_os << GENOTYPE_STATUS::label(obs.variantStatus) << "\t";
                obs_os << context.repeatCount << "\t" << INDEL_SIGNAL_TYPE::label(sigIndex) << "\t";
                obs_os << indelKey.deletionLength << "\t" << nonrefAlleleIndex + 1 << "/" << nonrefAlleleCount << "\t";
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

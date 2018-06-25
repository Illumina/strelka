//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2018 Illumina, Inc.
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

#include "SequenceAlleleCountsPosProcessor.hh"
#include "blt_common/ref_context.hh"
#include "common/OutStream.hh"
#include "common/Exceptions.hh"
#include "starling_common/OrthogonalVariantAlleleCandidateGroupUtil.hh"
#include "strelka_common/StrelkaSampleSetSummary.hh"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>



/// Checks that the path given in \p fileName can be opened for writing
///
/// If \p fileName is empty it is ignored, if \p fileName is not writable for any reason
/// (path does not exist, user does not have permission, etc..) this function will throw.
static
void
checkOutputFilePathIsWriteable(
    const std::string& fileName)
{
    OutStream outs(fileName);
}



SequenceAlleleCountsPosProcessor::
SequenceAlleleCountsPosProcessor(
    const SequenceAlleleCountsOptions& opt,
    const SequenceAlleleCountsDerivOptions& dopt,
    const reference_contig_segment& ref,
    const SequenceAlleleCountsStreams& fileStreams,
    RunStatsManager& statsManager)
    : base_t(opt, dopt, ref, fileStreams, opt.alignFileOpt.alignmentFilenames.size(), statsManager),
      _opt(opt),
      _dopt(dopt),
      _streams(fileStreams)
{
    // not generalized to multi-sample yet:
    assert(getSampleCount()==1);
    static const unsigned sampleId(0);

    // set the sample name to the first bam file name
    // this is already logged as a warning if alignmentFilenames.size() > 1
    _counts.setSampleName(opt.alignFileOpt.alignmentFilenames[0]);

    // check that we have write permission on all output files as early as possible:
    checkOutputFilePathIsWriteable(opt.countsFilename);
    if (opt.is_write_observations())
    {
        checkOutputFilePathIsWriteable(opt.observationsBedFilename);
    }
    if (! opt.nonEmptySiteCountFilename.empty())
    {
        checkOutputFilePathIsWriteable(opt.nonEmptySiteCountFilename);
    }

    // setup indel buffer samples
    {
        sample_info& normal_sif(sample(sampleId));

        const unsigned syncSampleId = getIndelBuffer().registerSample(normal_sif.estdepth_buff, normal_sif.estdepth_buff_tier2, true);

        assert(syncSampleId == sampleId);

        getIndelBuffer().finalizeSamples();
    }
}



void
SequenceAlleleCountsPosProcessor::
reset()
{
    base_t::reset();

    _excludedRegions.clear();
    _knownVariants.clear();
}



void
SequenceAlleleCountsPosProcessor::
completeProcessing()
{
    reset();

    _counts.save(_opt.countsFilename.c_str());

    if (! _opt.nonEmptySiteCountFilename.empty())
    {
        OutStream outs(_opt.nonEmptySiteCountFilename);
        outs.getStream() << "nonEmptySiteCount\t" << _nonEmptySiteCount << "\n";
    }

    _nonEmptySiteCount = 0;
}



void
SequenceAlleleCountsPosProcessor::
resetRegion(
    const std::string& chromName,
    const known_pos_range2& reportRegion)
{
    base_t::resetRegionBase(chromName, reportRegion);

    // setup norm and max filtration depths
    {
        if (_opt.is_depth_filter())
        {
            cdmap_t::const_iterator cdi(_dopt.chrom_depth.find(chromName));
            if (cdi == _dopt.chrom_depth.end())
            {
                std::ostringstream oss;
                oss << "Can't find chromosome: '" << chromName << "' in chrom depth file: '"
                    << _opt.chrom_depth_file << "'";
                throw blt_exception(oss.str().c_str());
            }
            _normChromDepth = cdi->second;
            _maxChromDepth = (_normChromDepth * _opt.max_depth_factor);
        }
        assert(_normChromDepth >= 0.);
        assert(_maxChromDepth >= 0.);
    }

    // setup indel buffer max depth
    {
        if (_dopt.is_max_depth())
        {
            if (_opt.max_candidate_indel_depth_factor > 0.)
            {
                _maxNormalSampleDepthForCandidateVariants = (_opt.max_candidate_indel_depth_factor * _maxChromDepth);
            }
        }

        if (_opt.max_candidate_indel_depth > 0.)
        {
            if (_maxNormalSampleDepthForCandidateVariants > 0.)
            {
                _maxNormalSampleDepthForCandidateVariants = std::min(_maxNormalSampleDepthForCandidateVariants,
                                                                     static_cast<double>(_opt.max_candidate_indel_depth));
            }
            else
            {
                _maxNormalSampleDepthForCandidateVariants = _opt.max_candidate_indel_depth;
            }
        }

        getIndelBuffer().setMaxCandidateDepth(_maxNormalSampleDepthForCandidateVariants);
    }
}



void
SequenceAlleleCountsPosProcessor::
insertExcludedRegion(
    const known_pos_range2& excludedRange)
{
    _stagemanPtr->validate_new_pos_value(excludedRange.begin_pos(),STAGE::READ_BUFFER);
    _excludedRegions.addRegion(excludedRange);
    _is_skip_process_pos=false;
}



void
SequenceAlleleCountsPosProcessor::
addKnownVariant(
    const vcf_record& knownVariant)
{
    _stagemanPtr->validate_new_pos_value(knownVariant.pos, STAGE::READ_BUFFER);
    _knownVariants.addVcfRecord(knownVariant);
    _is_skip_process_pos=false;
}


/// \brief test if the given indel key has a match in the knownVariants set
///
/// \param overlap populated with the matching known variant(s) if any are found
/// \return true if any known variant matches are found
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



/// return supporting read counts for all (overlapping) alleles in alleleGroup
///
/// \param support supporting read counts for the reference (on index 0), and for each allele in alleleGroup thereafter
static
void
getOrthogonalHaplotypeSupportCounts(
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    const unsigned sampleId,
    const unsigned minDistanceFromReadEdge,
    std::vector<unsigned>& support)
{
    static const bool isTier1Only(true);

    const unsigned nonrefAlleleCount(alleleGroup.size());
    assert(nonrefAlleleCount!=0);

    // intersection of read ids which have a likelihood evaluated over all candidate haplotypes:
    std::set<unsigned> readIds;
    getAlleleGroupIntersectionReadIds(sampleId, alleleGroup, readIds, isTier1Only, minDistanceFromReadEdge);

    // count of all haplotypes including reference
    const unsigned fullAlleleCount(nonrefAlleleCount+1);

    support.clear();
    support.resize(fullAlleleCount,0);

    for (const auto readId : readIds)
    {
        std::vector<double> lhood(fullAlleleCount);
        getAlleleNaivePosteriorFromRead(sampleId, alleleGroup, readId, lhood);
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
IndelCounts::INDEL_SIGNAL_TYPE::index_t
getIndelType(
    const AlleleReportInfo& indelReportInfo)
{
    int rudiff(static_cast<int>(indelReportInfo.refRepeatCount)-static_cast<int>(indelReportInfo.indelRepeatCount));
    rudiff = std::min(3,std::max(-3,rudiff));

    using namespace IndelCounts::INDEL_SIGNAL_TYPE;
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



/// Merge indel observations for a given position and context
///
/// Indel observations at the same position and context might need to be summed if these are different "signal"
/// types (where signal types are 1 base deletion, 2 base deletion, etc...)
///
static
void
mergeIndelObservations(
    const IndelCounts::Context& context,
    const IndelCounts::SingleSampleCandidateVariantContextObservationPattern& indelObservation,
    std::map<IndelCounts::Context,
    IndelCounts::SingleSampleCandidateVariantContextObservationPattern>& mergedIndelObservations)
{
    using namespace IndelCounts;

    auto iter(mergedIndelObservations.find(context));

    if (iter == mergedIndelObservations.end())
    {
        mergedIndelObservations.insert(std::make_pair(context,indelObservation));
    }
    else
    {
        // all signal counts should be summed:
        for (unsigned signalIndex(0); signalIndex<INDEL_SIGNAL_TYPE::SIZE; ++signalIndex)
        {
            iter->second.signalCounts[signalIndex] += indelObservation.signalCounts[signalIndex];
        }

        // refCount used for the group is the lowest submitted:
        iter->second.refCount = std::min(iter->second.refCount, indelObservation.refCount);

        // known variant status is the maximum status within the set (e.g. if UNKNOWN
        // and VARIANT, site is variant)
        iter->second.variantStatus = std::max(iter->second.variantStatus, indelObservation.variantStatus);
    }
}



/// all STR context information derived from a specific reference position
struct ReferenceSTRContext
{
    bool isBaseInSTR = false;
    bool isBaseLeftEndOfSTR = false;
    unsigned patternSize = 1;

    // not computed unless isBaseInSTR and isBaseLeftEndOfSTR
    unsigned STRRepeatCount = 1;
};



static const unsigned maxSTRRepeatCount(20);



/// Derive the STR context for a given position in the reference
static
ReferenceSTRContext
getReferenceSTRContext(
    const reference_contig_segment& ref,
    const pos_t pos)
{
    // STR pattern sizes to be counted
    static const std::vector<unsigned> referenceSTRPatternSizeVector = {1,2};

    ReferenceSTRContext refSTRContext;

    // find the pattern size of the STR track the current base is in
    // use the smaller pattern size if the base is in two STR tracks (e.g., pos 2 in AAAGAG is in the hpol track
    for (const auto patternSize : referenceSTRPatternSizeVector)
    {
        searchForSTR(patternSize, pos, refSTRContext.isBaseInSTR, refSTRContext.isBaseLeftEndOfSTR, ref);
        if (refSTRContext.isBaseInSTR)
        {
            refSTRContext.patternSize = patternSize;
            if (refSTRContext.isBaseLeftEndOfSTR)
            {
                refSTRContext.STRRepeatCount =
                    std::min(maxSTRRepeatCount, getLeftShiftedSTRRepeatCount(patternSize, pos, ref));
            }
            break;
        }
    }

    return refSTRContext;
}



void
SequenceAlleleCountsPosProcessor::
process_pos_error_counts(
    const pos_t pos)
{
    if (! is_pos_reportable(pos)) return;

    const unsigned sampleCount(getSampleCount());

    // the error counts workflow can only be called for a single sample:
    //
    assert(sampleCount==1);
    const unsigned sampleIndex(0);
    sample_info& sif(sample(sampleIndex));


    const char refBase(_ref.get_base(pos));

    if (refBase=='N') return;

    BasecallCounts::Dataset& basecallCounts(_counts.getBasecallCounts());
    IndelCounts::Dataset& indelCounts(_counts.getIndelCounts());

    // right now there's only one basecallContext (ie. only one context used for SNVs and basecalls), this
    // is why we set it as const here
    //
    // this contrasts with the indel counting logic below which is using several contexts
    const BasecallCounts::Context basecallContext;

    bool isSkipSNV(false);
    bool isSkipIndel(false);

    const ReferenceSTRContext referenceSTRContext(getReferenceSTRContext(_ref,pos));

    // For consistency, we don't take any indel counts information from an STR track except for the left-most position
    //
    // Note that we may want to change this in the future because calls can occur within STRs which are not expansions
    // or contractions of the STR pattern, and thus will not be left-shifted, e.g. "AAAA -> AACAA"
    //
    isSkipIndel = (referenceSTRContext.isBaseInSTR and (not referenceSTRContext.isBaseLeftEndOfSTR));


    if (_excludedRegions.isIntersectRegion(pos))
    {
        // The workflow has the option to exclude some regions, so count the number of excluded sites.
        // Excluded sites are recorded for QC purposes but don't impact rate estimates.
        basecallCounts.addExcludedRegionSkip(basecallContext);
        isSkipSNV=true;

        // record the excluded region count only for positions where we would have counted indel evidence anyway:
        if (not isSkipIndel)
        {
            IndelCounts::Context indelContext(referenceSTRContext.patternSize, referenceSTRContext.STRRepeatCount);
            indelCounts.addExcludedRegionSkip(indelContext);
            isSkipIndel=true;
        }
    }
    else if (_maxNormalSampleDepthForCandidateVariants > 0.)
    {
        // determine if SNV or indel output is going to be disabled because we're in a region of anomalously high depth
        //
        // Candidate variants are turned off above a certain depth relative to chromosome mean.
        // For basecalls we check the relative depth at current position, for indels we check the relative depth at previous position.
        //
        {
            // handle SNVs
            const unsigned estdepth(sif.estdepth_buff.val(pos));
            const unsigned estdepth2(sif.estdepth_buff_tier2.val(pos));
            if ((estdepth+estdepth2) > _maxNormalSampleDepthForCandidateVariants)
            {
                basecallCounts.addDepthSkip(basecallContext);
                isSkipSNV=true;
            }
        }

        {
            const unsigned estdepth(sif.estdepth_buff.val(pos-1));
            const unsigned estdepth2(sif.estdepth_buff_tier2.val(pos-1));
            if ((estdepth+estdepth2) > _maxNormalSampleDepthForCandidateVariants)
            {
                // record the high-depth bypass count only for positions where we would have counted indel evidence anyway:
                if (not isSkipIndel)
                {
                    IndelCounts::Context indelContext(referenceSTRContext.patternSize, referenceSTRContext.STRRepeatCount);
                    indelCounts.addDepthSkip(indelContext);
                    isSkipIndel=true;
                }
            }
        }
    }


    // handle basecall error signal
    if (! isSkipSNV)
    {
        // even though this site may still be excluded from SNV counting because of noise, we use this counter to
        // approximate how much of the genome has some level of detectable coverage that could possibly contribute to
        // error estimation
        _nonEmptySiteCount += 1;

        // be relatively intolerant of anything interesting happening in the local sequence neighborhood:
        static const double snvMMDRMaxFrac(0.05);

        const snp_pos_info& sinfo(sif.basecallBuffer.get_pos(pos));
        unsigned fcount(0);
        for (const base_call& bc : sinfo.calls)
        {
            if (bc.is_call_filter) fcount++;
        }
        fcount += sinfo.tier2_calls.size();

        const double filtFrac(fcount/static_cast<double>(sinfo.calls.size()+sinfo.tier2_calls.size()));

        if (filtFrac >= snvMMDRMaxFrac)
        {
            // skip basecall stats at this position because of a noisy local environment
            //
            // track count of how often this happens:
            basecallCounts.addNoiseSkip(basecallContext);
        }
        else
        {
            bool isEmpty(true);
            const uint8_t ref_id(base_to_id(refBase));
            BasecallCounts::ContextInstanceObservation observation;
            for (const base_call& bc : sinfo.calls)
            {
                if (bc.is_call_filter) continue;

                const uint16_t qual(bc.get_qscore());
                if (qual<_opt.minBasecallErrorPhredProb) continue;

                static const uint16_t min_count_qscore(25);
                if (qual<min_count_qscore) continue;
                isEmpty=false;
                const bool isRef(bc.base_id == ref_id);
                if (isRef)
                {
                    observation.addRefCount(bc.is_fwd_strand, bc.get_qscore());
                }
                else
                {
                    observation.addAltCount(bc.is_fwd_strand, bc.get_qscore());
                }
            }

            if (isEmpty)
            {
                basecallCounts.addEmptySkip(basecallContext);
            }
            else
            {
                basecallCounts.addContextInstanceObservation(basecallContext, observation);
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
    // large. Otherwise, read support for each allele in the overlapping set is enumerated. For all alleles that
    // start at this position, we record the supporting counts in the indel error analysis counting structure.
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
            bool isZeroCoverage(true);
            for (unsigned sampleIndex2(0); sampleIndex2 < sampleCount; ++sampleIndex2)
            {
                const IndelSampleData& indelSampleData(indelData.getSampleData(sampleIndex2));
                if (not indelSampleData.read_path_lnp.empty())
                {
                    isZeroCoverage = false;
                    break;
                }
            }

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
    // in this case we choose the latter in principal (favor true noise vs what the germline caller sees),
    // but still set an upper-limit on the total number of overlapping variants, recognizing that we can't
    // handle the noise accurately as variant density goes up beyond a certain point.
    //
    // Future design notes:
    // This parameter might have to be revisited for long homopolymers in particular, where both mutation and
    // sequencing error rates are high enough that greater than 4 alt alleles may be a reasonably common occurance.
    //
    // Alternatively -- we could take the top maxOverlap alleles instead of skipping the locus entirely, as we do now
    //
    const unsigned maxOverlap(4);
    std::vector<unsigned> callerPloidy = { maxOverlap };


    if (orthogonalVariantAlleles.size() > maxOverlap) return;

    // check for any known variants overlapping this position
    RecordTracker::indel_value_t knownVariantRecords;
    _knownVariants.intersectingRecord(pos, knownVariantRecords);

    // buffer observations (as SingleSampleContextObservationPattern objects) until we get through all overlapping indels
    // at this position, then process these into the error counts accumulation data structure:
    //
    // this buffering allows us to cleanly break out of this locus if we see evidence of too much noise below
    //
    std::map<IndelCounts::Context,
        IndelCounts::SingleSampleCandidateVariantContextObservationPattern> mergedIndelObservations;


    if (not orthogonalVariantAlleles.empty())
    {
        {
            std::vector<unsigned> topVariantAlleleIndexPerSample;
            const bool isEveryAltIncluded = addAllelesAtOtherPositions(
                                                _ref, sampleCount, callerPloidy, pos, get_largest_total_indel_ref_span_per_read(),
                                                getIndelBuffer(), orthogonalVariantAlleles, topVariantAlleleIndexPerSample);

            if (not isEveryAltIncluded) return;
        }

        if (orthogonalVariantAlleles.size() > maxOverlap) return;

        const unsigned nonrefAlleleCount(orthogonalVariantAlleles.size());
        std::vector<unsigned> support;
        getOrthogonalHaplotypeSupportCounts(
            orthogonalVariantAlleles, sampleIndex, _opt.minDistanceFromReadEdge, support);

        for (unsigned nonrefAlleleIndex(0); nonrefAlleleIndex < nonrefAlleleCount; ++nonrefAlleleIndex)
        {
            const IndelKey& indelKey(orthogonalVariantAlleles.key(nonrefAlleleIndex));

            if (indelKey.pos != pos) continue;

            const IndelData& indelData(orthogonalVariantAlleles.data(nonrefAlleleIndex));
            const AlleleReportInfo& indelReportInfo(indelData.getReportInfo());

            IndelCounts::Context context;

            if ((indelReportInfo.repeatUnitLength == referenceSTRContext.patternSize) && (indelReportInfo.refRepeatCount > 1))
            {
                // guard against the occasional non-normalized indel:
                if (referenceSTRContext.STRRepeatCount == std::min(maxSTRRepeatCount, indelReportInfo.refRepeatCount))
                {
                    context = IndelCounts::Context(
                                  referenceSTRContext.patternSize, referenceSTRContext.STRRepeatCount);
                }
            }


            // check to see if this indel is (likely to be) a match to a known variant
            RecordTracker::indel_value_t overlappingRecords;
            isKnownVariantMatch(knownVariantRecords, indelKey, overlappingRecords);

            IndelCounts::SingleSampleCandidateVariantContextObservationPattern indelObservation;

            const IndelCounts::INDEL_SIGNAL_TYPE::index_t sigIndex(getIndelType(indelReportInfo));
            indelObservation.signalCounts[sigIndex] = support[nonrefAlleleIndex + 1];
            static const unsigned refAlleleIndex(0);
            indelObservation.refCount = support[refAlleleIndex];
            indelObservation.assignKnownStatus(overlappingRecords);

            // debug output:
            //
            // an indel candidate can have 0 q30 indel reads when it is only supported by
            // noise reads (i.e. indel occurs outside of a read's valid alignment range,
            // see lib/starling_common/starling_read_util.cpp::get_valid_alignment_range)
            // in this case, we're not going to report the incidence as noise, since it's
            // not a read we would consider in variant calling
            if (support[nonrefAlleleIndex + 1] > 0 && _opt.is_write_observations())
            {
                std::ostream& obs_os(*_streams.observation_bed_osptr());
                obs_os << _chromName << "\t";
                obs_os << indelKey.pos << "\t" << indelKey.pos + indelKey.deletionLength << "\t"
                       << INDEL::get_index_label(indelKey.type) << "\t";
                obs_os << indelReportInfo.repeatUnit << "\t" << indelReportInfo.refRepeatCount << "\t";
                obs_os << GENOTYPE_STATUS::label(indelObservation.variantStatus) << "\t";
                obs_os << context.getRepeatPatternSize() << "\t" << context.getRepeatCount() << "\t"
                       << IndelCounts::INDEL_SIGNAL_TYPE::label(sigIndex) << "\t";
                obs_os << indelKey.deletionLength << "\t" << nonrefAlleleIndex + 1 << "/" << nonrefAlleleCount
                       << "\t";
                obs_os << indelObservation.signalCounts[sigIndex] << "\t" << indelObservation.refCount << "\t";
                obs_os << std::accumulate(support.begin(), support.end(), 0) << std::endl;

#if 0
                for (const auto& rec : overlappingRecords)
                {
                    std::cout << rec << std::endl;
                }
#endif
            }
            mergeIndelObservations(context, indelObservation, mergedIndelObservations);
        }
    }



    // background depth is always one minus position to be consistent with indel report:
    const pos_t depth_pos(pos - 1);
    const snp_pos_info& spi(sif.basecallBuffer.get_pos(depth_pos));
    const unsigned depth(spi.calls.size());

    for (const auto& value : mergedIndelObservations)
    {
        indelCounts.addCandidateVariantContextInstanceObservation(value.first, value.second, depth);
    }

    // add all the contexts that haven't been covered already by indel candidate variants
    {
        IndelCounts::Context context;

        if (referenceSTRContext.isBaseInSTR)
        {
            context = IndelCounts::Context(referenceSTRContext.patternSize, referenceSTRContext.STRRepeatCount);
        }

        IndelCounts::SingleSampleNonVariantContextObservationPattern nonVariantObservation;
        nonVariantObservation.depth = depth;

        // the assumption is that a background position should have
        // the variant status of any overlapping known variants,
        // regardless of whether the genotypes match
        nonVariantObservation.assignKnownStatus(knownVariantRecords);

        if (! mergedIndelObservations.count(context))
        {
            indelCounts.addNonVariantContextInstanceObservation(context, nonVariantObservation);
        }
    }
}

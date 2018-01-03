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

///
/// \author Chris Saunders
///

#include "alignment_util.hh"
#include "starling_read_align_score_indels.hh"

#include "blt_util/log.hh"
#include "starling_common/indel_util.hh"

#include <algorithm>
#include <iostream>


//#define DEBUG_ALIGN

// iks_map_t stores information linking candidate indels to alignment scores
//
// The components of the map type are enumerated below. The "active" indel
// is the one being evaluated/scored.
//

/// first: is X indel present in the read alignment?
/// second: X indel, this may or may not be the same as the active indel,
///         but it should only be the active indel or an indel which interferes
///         with the active indel
typedef std::pair<bool,IndelKey> indel_present_t;

/// first: active indel
typedef std::pair<IndelKey, indel_present_t > indel_status_t;

/// first: alignment likelihood
/// second: alignment
typedef std::pair<double,const CandidateAlignment*> align_info_t;
typedef std::map<indel_status_t,align_info_t> iks_map_t;



/// check whether the current alignment has a better score than that already
/// stored in iks_map and if so add entry to iks_map
///
static
void
updateIndelScoringInfo(
    iks_map_t& iks_map,
    const IndelKey& ik_call,
    const bool is_indel_present,
    const IndelKey& ik_present,
    const double path_lnp,
    const CandidateAlignment* cal_ptr)
{
    const indel_status_t mkey(std::make_pair(ik_call,std::make_pair(is_indel_present,ik_present)));

    // check to see if the path score is better than what we already have:
    const iks_map_t::const_iterator j(iks_map.find(mkey));
    if ((j!=iks_map.end()) && ((j->second.first) >= path_lnp)) return;
    iks_map[mkey]=std::make_pair(path_lnp,cal_ptr);
}



typedef std::map<IndelKey,indel_set_t> overlap_map_t;



static
void
overlap_map_tick(
    overlap_map_t& omap,
    const IndelKey& indelKey1,
    const IndelKey& indelKey2)
{
    indel_set_t& os(omap[indelKey1]);
    if (os.find(indelKey2) == os.end()) os.insert(indelKey2);
}



/// if "new_indel" interferes with any of the indel set "current_indels", return iterator pointing to one
/// current_indels member which is conflicing
static
indel_set_t::const_iterator
which_interfering_indel(
    const indel_set_t& current_indels,
    const IndelKey& new_indel)
{
    typedef indel_set_t::const_iterator ci;
    ci end(current_indels.end());
    if (current_indels.count(new_indel) != 0) return end;

    ci iter(current_indels.begin());
    const ci iter_end(current_indels.end());
    for (; iter != iter_end; ++iter)
    {
        const IndelKey& currentIndel(*iter);
        if (currentIndel.isMismatch()) continue;
        if (is_indel_conflict(currentIndel, new_indel)) return iter;
    }
    return end;
}



/// by how many bases does the alignment cross each indel breakpoint?
///
/// Note that this function only works correctly for indels that are
/// present in the alignment -- otherwise there are multiple possible
/// realignments of the read
///
/// Allow for the possibility of an upstream oligo to be implicitly
/// present for the read, acting to 'anchor' the beginning of the read
/// with the equivalent length match state even though this does not
/// appear in the read's CIGAR string.
///
static
std::pair<int,int>
get_alignment_indel_bp_overlap(
    const unsigned upstream_oligo_size,
    const alignment& al,
    const IndelKey& indelKey)
{
    using namespace ALIGNPATH;

    // all read positions are expressed in unclipped read coordinates!
    //
    pos_t read_head_pos(0);
    pos_t ref_head_pos(al.pos);

    // mark the left and right indel breakpoints in read position
    // coordinates:
    bool is_left_read_pos(false);
    pos_t left_read_pos(0);
    bool is_right_read_pos(false);
    pos_t right_read_pos(0);

    for (const path_segment& ps : al.path)
    {
        pos_t next_read_head_pos(read_head_pos);
        pos_t next_ref_head_pos(ref_head_pos);

        if       (is_segment_align_match(ps.type))
        {
            next_read_head_pos += ps.length;
            next_ref_head_pos += ps.length;
        }
        else if (ps.type == INSERT)
        {
            next_read_head_pos += ps.length;
        }
        else if (ps.type == DELETE)
        {
            next_ref_head_pos += ps.length;
        }
        else if ((ps.type == SOFT_CLIP) || (ps.type == HARD_CLIP))
        {
            // do nothing... this function operates on unclipped read
            // coordinates instead of true read coordinates
            //
        }
        else
        {
            // TODO deal with other CIGAR types....
            assert(0);
        }

        if ((! is_left_read_pos) && (indelKey.pos<=(next_ref_head_pos)))
        {
            left_read_pos=read_head_pos+(indelKey.pos-ref_head_pos);
            is_left_read_pos=true;
        }
        if ((! is_right_read_pos) && (indelKey.right_pos()<(next_ref_head_pos)))
        {
            right_read_pos=read_head_pos+(indelKey.right_pos()-ref_head_pos);
            is_right_read_pos=true;
        }

        read_head_pos=next_read_head_pos;
        ref_head_pos=next_ref_head_pos;
    }

    // extension values account for oligo anchors at the end of the
    // read not represented in the CIGAR string:
    int left_extension(0);
    int right_extension(0);
    if (al.is_fwd_strand)
    {
        if ( left_read_pos > 0 ) left_extension = upstream_oligo_size;
    }
    else
    {
        if ( (read_head_pos-right_read_pos)  > 0 )right_extension = upstream_oligo_size;
    }

    int left_overlap(0);
    if (is_left_read_pos)
    {
        left_overlap=std::min(left_read_pos+left_extension,(read_head_pos-left_read_pos));
        left_overlap=std::max(0,left_overlap);
    }

    int right_overlap(0);
    if (is_right_read_pos)
    {
        right_overlap=std::min(right_read_pos,(read_head_pos-right_read_pos)+right_extension);
        right_overlap=std::max(0,right_overlap);
    }

    return std::make_pair(left_overlap,right_overlap);
}



typedef std::set<std::pair<IndelKey,IndelKey> > indel_pair_set;



/// given two alignments that have (nearly) the same alignment score,
/// quickly approximate whether they could be expressing the same set
/// of indels:
///
static
bool
is_equiv_candidate(
    const CandidateAlignment& cal1,
    const CandidateAlignment& cal2,
    indel_pair_set& equiv_keys)
{
    equiv_keys.clear();

    const indel_set_t& is1(cal1.getIndels());
    const indel_set_t& is2(cal2.getIndels());

    const unsigned s1(is1.size());
    const unsigned s2(is2.size());

    if (s1 != s2) return false;

    indel_set_t::const_iterator it1(is1.begin()), it1_end(is1.end());
    indel_set_t::const_iterator it2(is2.begin()); //, it2_end(is2.end());

    for (; it1!=it1_end; ++it1,++it2)
    {
        if (*it1==*it2) continue;
        if (it1->type != it2->type) return false;
        if (it1->delete_length() != it2->delete_length()) return false;
        if (it1->insert_seq() != it2->insert_seq()) return false;
        equiv_keys.insert(std::make_pair(*it1,*it2));
    }
    return true;
}



/// written for equiv indel resolution (i.e. assuming ik1 and ik2 are
/// of the same type and length)
///
static
bool
is_first_indel_dominant(
    const IndelBuffer& indelBuffer,
    const IndelKey& indelKey1,
    const IndelKey& indelKey2)
{
    const bool ic1(indelBuffer.isCandidateIndel(indelKey1));
    const bool ic2(indelBuffer.isCandidateIndel(indelKey2));

    if (ic2 && (! ic1)) return false;
    if (ic2==ic1)
    {
        return (indelKey1.pos<=indelKey2.pos);
    }
    return true;
}



/// A heuristic late-stage indel normalization procedure,
/// which identifies indels likely to represent the same
/// haplotype and disables all but the first.
///
/// \param[out] nonnorm_indels disabled indel set
/// \param[out] isFilterCandAlignment true if alignment contains an excluded indel
///
static
void
late_indel_normalization_filter(
    const starling_base_options& opt,
    const IndelBuffer& indelBuffer,
    const std::set<CandidateAlignment>& candAlignments,
    const std::vector<double>& candAlignmentScores,
    indel_set_t nonnorm_indels,
    std::vector<bool>& isFilterCandAlignment,
    double& maxCandAlignmentScore,
    const CandidateAlignment*& maxCandAlignmentPtr)
{
    const unsigned candAlignmentCount(candAlignments.size());

    const double equiv_lnp_range( opt.is_smoothed_alignments ?
                                  opt.smoothed_lnp_range : 0. );

    // go through alignment x alignments in best->worst score
    // order -- to do so start out with the sort order:
    std::vector<std::pair<double,unsigned> > sortedScores;
    std::vector<const CandidateAlignment*> candAlignmentPtrs;
    {
        std::set<CandidateAlignment>::const_iterator si(candAlignments.begin());
        for (unsigned candAlignmentIndex(0); candAlignmentIndex<candAlignmentCount; ++candAlignmentIndex,++si)
        {
            sortedScores.push_back(std::make_pair(candAlignmentScores[candAlignmentIndex],candAlignmentIndex));
            candAlignmentPtrs.push_back(&(*si));
        }
        std::sort(sortedScores.rbegin(),sortedScores.rend());
    }

    // this starts out as a copy of candAlignmentScores, but smoothScores
    // gets reset to the highest lnp an alignment has eliminated
    // to prevent transitive shenanigans:
    std::vector<double> smoothScores(candAlignmentScores);

    bool is_any_excluded(false);
    indel_pair_set indelPairs;

    for (unsigned index1(0); index1<candAlignmentCount; ++index1)
    {
        const unsigned sortedIndex1(sortedScores[index1].second);
        if (isFilterCandAlignment[sortedIndex1]) continue;

        for (unsigned index2(index1+1); index2<candAlignmentCount; ++index2)
        {
            const unsigned sortedIndex2(sortedScores[index2].second);
            if (isFilterCandAlignment[sortedIndex2]) continue;

            // determine if this is an equiv_indel:
            if (smoothScores[sortedIndex2]+equiv_lnp_range < smoothScores[sortedIndex1]) break;
            const bool is_equiv(is_equiv_candidate(*(candAlignmentPtrs[sortedIndex1]),
                                                   *(candAlignmentPtrs[sortedIndex2]),
                                                   indelPairs));

            if (! is_equiv) continue;

#ifdef DEBUG_ALIGN
            log_os << "COWSLIP: sorti,lnpi,sortj,lnpj: "
                   << sortedIndex1 << " " << smoothScores[sortedIndex1] << " "
                   << sortedIndex2 << " " << smoothScores[sortedIndex2] << "\n";
            log_os << "COWSLIP cali: " << *(candAlignmentPtrs[sortedIndex1]) << "\n";
            log_os << "COWSLIP calj: " << *(candAlignmentPtrs[sortedIndex2]) << "\n";
#endif

            // In theory, this should be an assertion of
            // !ips.empty(), but it seems that two alignments with
            // the same score and the same set of indels, but
            // (i.e. different start positions) could somehow
            // occur.
            if (indelPairs.empty()) continue;

            // for each pair, determine which indel is dominant,
            // if there are multiple equiv indels, we remove the
            // alignment containing the first non-dominant equiv
            // indel only
            bool is_sortedIndex1_removed(false);
            bool is_removed(false);
            for (const auto& indelPair : indelPairs)
            {
                const bool is1(is_first_indel_dominant(indelBuffer, indelPair.first, indelPair.second));

#ifdef DEBUG_ALIGN
                log_os << "COWSLIP: indel1: " << indelPair.first << "\n";
                log_os << "COWSLIP: indel2: " << indelPair.second << "\n";
                log_os << "COWSLIP: is_indel1_dominant?: " << is1 << "\n";
#endif
                if (is1)
                {
                    nonnorm_indels.insert(indelPair.second);
#ifdef DEBUG_ALIGN
                    log_os << "COWSLIP: marking 2 nonnorm: " << indelPair.second << "\n";
#endif
                    if (! is_removed)
                    {
                        isFilterCandAlignment[sortedIndex2] = true;
                        is_any_excluded = true;
                        smoothScores[sortedIndex1] = std::max(smoothScores[sortedIndex1],smoothScores[sortedIndex2]);
#ifdef DEBUG_ALIGN
                        log_os << "COWSLIP: excluding sortj: " << sortedIndex2 << "\n";
#endif
                    }
                }
                else
                {
                    nonnorm_indels.insert(indelPair.first);
#ifdef DEBUG_ALIGN
                    log_os << "COWSLIP: marking 1 nonnorm: " << indelPair.first << "\n";
#endif
                    if (! is_removed)
                    {
                        isFilterCandAlignment[sortedIndex1] = true;
                        is_any_excluded = true;
                        smoothScores[sortedIndex2] = std::max(smoothScores[sortedIndex1],smoothScores[sortedIndex2]);
#ifdef DEBUG_ALIGN
                        log_os << "COWSLIP: excluding sorti: " << sortedIndex1 << "\n";
#endif
                        is_sortedIndex1_removed=true;
                    }
                }
                is_removed=true;
            }
            if (is_sortedIndex1_removed) break;
        }
    }

    // recalc new max path in case it was excluded:
    //
    // note this loop is designed to take advantage of the high->low path_lnp sort
    //
    if (is_any_excluded)
    {
        for (unsigned candAlignmentIndex(0); candAlignmentIndex<candAlignmentCount; ++candAlignmentIndex)
        {
            const unsigned sortedIndex(sortedScores[candAlignmentIndex].second);
            if (isFilterCandAlignment[sortedIndex]) continue;
#ifdef DEBUG_ALIGN
            log_os << "COWSLIP: reseting max_path_lnp/max_cal from: " << maxCandAlignmentScore << " " << *(maxCandAlignmentPtr) << "\n";
#endif
            maxCandAlignmentScore=candAlignmentScores[sortedIndex];
            maxCandAlignmentPtr=candAlignmentPtrs[sortedIndex];
#ifdef DEBUG_ALIGN
            log_os << "COWSLIP: reseting max_path_lnp/max_cal to: " << maxCandAlignmentScore << " " << *(maxCandAlignmentPtr) << "\n";
#endif
            break;
        }
    }
}



void
score_indels(
    const starling_base_options& opt,
    const starling_base_deriv_options&,
    const starling_sample_options& sample_opt,
    const read_segment& rseg,
    IndelBuffer& indelBuffer,
    const unsigned sampleIndex,
    const std::set<CandidateAlignment>& candAlignments,
    const bool is_incomplete_search,
    const std::vector<double>& candAlignmentScores,
    double maxCandAlignmentScore,
    const CandidateAlignment* maxCandAlignmentPtr)
{
    static const bool is_safe_mode(true);

    // (1) score candidate alignments -- already done before calling this function
    //
    // (2) find the highest scoring alignment with each indel present
    // and absent.
    //
    // Note there are up to two 'indel absent' states. The first is
    // when an interfering indel is not present, and the second is
    // when an interfering indel is present. The second case will be
    // less common.
    //
    // We must control for the case where the indel is absent but
    // replaced by an equivalent indel. In this case we must not only
    // determine that an equivalent exists, but establish which case
    // is dominant so that the indel is only reported once. Dominance
    // rules should follow the same pattern established for
    // realignment above.
    //
    // The indel normalization here is intended to be a final fallback
    // corrective measure to prevent indels from disappearing at the
    // call stage. It is expected (and preferred) that as much
    // upstream indel normalization will take place as is practical.
    //

    // (2alpha)
    //
    // Go through all scored alignments and search for cases with an
    // equal score (or within the smooth threshold if this is turned
    // on). Determine if smooth score matches are actually equivalent
    // -- in which case remove the non-dominant alignment and mark the
    // equivalent indel as not participating in this read's indel call
    // search:
    //
    const unsigned candAlignmentCount(candAlignments.size());

    // Nonnorm indels contains any indels found to be equivalent to
    // another indel (ie. not normalized). These indels are not used
    // in subsequent calculation so as not to penalize the normalized
    // indel form (otherwise all reads which support the indel would
    // appear to ambiguously support two different indels).
    //
    indel_set_t nonnorm_indels;
    std::vector<bool> isFilterCandAlignment(candAlignmentCount,false);

    // This procedure implements a heuristic late-stage indel
    // normalization, set is_slip_norm=false to disable:
    //
    static const bool is_slip_norm(true);
    if (is_slip_norm)
    {
        late_indel_normalization_filter(opt, indelBuffer, candAlignments, candAlignmentScores, nonnorm_indels,
                                        isFilterCandAlignment, maxCandAlignmentScore, maxCandAlignmentPtr);
    }

    // (2a) find set of candidate indels which will be scored by this read
    //
    // The criteria for indels in this set are:
    // 1) the indel is not private to the read
    // 2) the indel is not in the nonnorm set
    // 3) the read overlaps at least one indel breakpoint by at least option:min_read_bp_flank bases in at least one candidate alignment
    //

    const CandidateAlignment& maxCandAlignment(*maxCandAlignmentPtr);

    indel_set_t indelsToEvaluate;
    const indel_set_t& indelsInMaxCandAlignment(maxCandAlignment.getIndels());
    {
        const known_pos_range maxCandAlignmentRange(get_soft_clip_alignment_range(maxCandAlignment.al));

#ifdef DEBUG_ALIGN
        log_os << "VARMIT: starting max_cal_eval_indel search\n"
               << "VARMIT: max_cal.al: " << maxCandAlignment.al << "\n"
               << "VARMIT: max_pr: " << maxCandAlignmentRange << "\n";
#endif

#ifdef DEBUG_ALIGN
        log_os << "VARMIT max_path extracted indels:\n";
        dump_indel_set(indelsInMaxCandAlignment,log_os);
#endif
        const auto indelIterPair(indelBuffer.rangeIterator(maxCandAlignmentRange.begin_pos, maxCandAlignmentRange.end_pos));
        for (auto indelIter(indelIterPair.first); indelIter!=indelIterPair.second; ++indelIter)
        {
            const IndelKey& evaluationIndel(indelIter->first);

            // mismatches are not scored
            if (evaluationIndel.isMismatch()) continue;

            const IndelData& indelData(getIndelData(indelIter));

#ifdef DEBUG_ALIGN
            log_os << "VARMIT: max path eval indel candidate: " << evaluationIndel;
#endif

            if (not indelBuffer.isCandidateIndel(evaluationIndel, indelData)) continue;
            if (nonnorm_indels.count(evaluationIndel) != 0) continue;

#ifdef DEBUG_ALIGN
            log_os << "VARMIT: max path indel is candidate and normed\n";
#endif

            const bool isIndelInMaxCandAlignment(indelsInMaxCandAlignment.count(evaluationIndel)!=0);
#ifdef DEBUG_ALIGN
            log_os << "VARMIT: indel present? " << isIndelInMaxCandAlignment << "\n";
#endif

            // (1) find the most likely alignment which contains this indel
            // (2) on that alignment, test for breakpoint overlap
            // (3) filter the candidate indel if overlap < min
            {
                const CandidateAlignment* maxCandAlignmentForIndelPtr(nullptr);
                if (isIndelInMaxCandAlignment)
                {
                    maxCandAlignmentForIndelPtr = &maxCandAlignment;
                }
                else
                {
                    double maxScore(0);
                    std::set<CandidateAlignment>::const_iterator candAlignmentIter(candAlignments.begin()),candAlignmentIter_end(candAlignments.end());
                    for (unsigned candAlignmentIndex(0); candAlignmentIter!=candAlignmentIter_end; ++candAlignmentIter,++candAlignmentIndex)
                    {
                        const CandidateAlignment& candAlignment(*candAlignmentIter);
                        const bool isMaxCandAlignment(&candAlignment == &maxCandAlignment);
                        if (isMaxCandAlignment) continue;

                        if (isFilterCandAlignment[candAlignmentIndex]) continue;

                        const indel_set_t& indelsInCandAlignment(candAlignment.getIndels());
                        const bool isIndelInCandAlignment(indelsInCandAlignment.count(evaluationIndel)!=0);
                        if (! isIndelInCandAlignment) continue;

                        const double score(candAlignmentScores[candAlignmentIndex]);
                        if ((nullptr == maxCandAlignmentForIndelPtr) || (score > maxScore))
                        {
                            maxScore = score;
                            maxCandAlignmentForIndelPtr = &(candAlignment);
                        }
                    }
                }

                // test and reject if breakpoint overlap is insufficient:
                if (nullptr == maxCandAlignmentForIndelPtr) continue;
                const std::pair<int,int> both_bpo(get_alignment_indel_bp_overlap(opt.upstream_oligo_size,maxCandAlignmentForIndelPtr->al,evaluationIndel));

                const int bpo(std::max(both_bpo.first,both_bpo.second));
#ifdef DEBUG_ALIGN
                log_os << "VARMIT: indel bp_overlap " << bpo << "\n";
#endif
                if (bpo < sample_opt.min_read_bp_flank)
                {
                    if (bpo>0)
                    {
                        const bool is_tier1_read(rseg.is_tier1_mapping());
                        IndelSampleData& indelSampleData(getIndelData(indelIter).getSampleData(sampleIndex));

                        if (is_tier1_read) indelSampleData.suboverlap_tier1_read_ids.insert(rseg.getReadIndex());
                        else               indelSampleData.suboverlap_tier2_read_ids.insert(rseg.getReadIndex());
                    }
                    continue;
                }
            }

            // all checks passed... indel will be evaluated for indel calling:
            //
#ifdef DEBUG_ALIGN
            log_os << "VARMIT: inserting indel to evaluate " << evaluationIndel << "\n";
#endif
            indelsToEvaluate.insert(evaluationIndel);
        }
    }

    // go through indelsToEvaluate and map which indels conflict with
    // each other
    //
    // for now we calc and store this info in a comically inefficient
    // manner
    //
    overlap_map_t orthogonalIndelMap;
    {
        indel_set_t::const_iterator i(indelsToEvaluate.begin());
        const indel_set_t::const_iterator i_end(indelsToEvaluate.end());
        for (; i != i_end; ++i)
        {
            indel_set_t::const_iterator j(i);
            ++j;
            for (; j != i_end; ++j)
            {
                if (is_indel_conflict(*i, *j))
                {
                    overlap_map_tick(orthogonalIndelMap, *i, *j);
                    overlap_map_tick(orthogonalIndelMap, *j, *i);
                }
            }
        }
    }

    // (2b) Create index of highest scoring alignment set for each
    //      evaluation indel (found in step 2a). The set of alignments
    //      include all mutually consistent sets of indels which
    //      interfere with the evaluation indel (including the eval
    //      indel itself) 99% of the time this should just be the
    //      indel and ref paths.  99% of the rest of the time this
    //      will just be the paths of the indel, ref, and one
    //      alternate indel allele.
    //

    // First part of iks_map_t key is the called indel -- second part
    // of the key is the indel which is toggled on which interferes
    // with the called indel, and a bool indicating the toggle state.
    // When called indel is toggled on the first and second indel keys
    // are the same. When bool is toggled off it is the no-indel
    // (reference) allele case.
    //
    // Note that we store maximum scoring alignment info for called indel,
    // reference, and each interfering indel allele. Most of the time there
    // should not be an interfering indel.
    //
    iks_map_t indelScoringInfo;
    {
        std::set<CandidateAlignment>::const_iterator candAlignmentIter(candAlignments.begin()),candAlignmentIter_end(candAlignments.end());
        for (unsigned candAlignmentIndex(0); candAlignmentIter!=candAlignmentIter_end; ++candAlignmentIter,++candAlignmentIndex)
        {
            const CandidateAlignment& candAlignment(*candAlignmentIter);
            const bool isMaxCandAlignment(&candAlignment == &maxCandAlignment);

            if (isFilterCandAlignment[candAlignmentIndex])
            {
                assert(! isMaxCandAlignment);
                continue;
            }

            const double score(candAlignmentScores[candAlignmentIndex]);

            const indel_set_t& indelsInCandAlignment(candAlignment.getIndels());

            indel_set_t nonCandidateIndelsOrthogonalToEvaluationIndels;

            //
            // SPECIAL NOTATION NOTE:
            //
            // In the following section there is a lot of documentation about a read having X haplotype in the sample,
            // but observed as Y haplotype, where the (X->Y) transformation represents a sequencing error.
            //
            // For consistency this is noted as:
            //
            // [[SAMPLE haplotype X TO OBSERVED haplotype Y]]
            //
            // This is intentionally overly-verbose b/c there have been past bugs and confusion on how we handle
            // these relationships
            //

            for (const IndelKey& evaluationIndel : indelsToEvaluate)
            {
                const IndelData* indelDataPtr(indelBuffer.getIndelDataPtr(evaluationIndel));
                assert(indelDataPtr != nullptr);
                const IndelSampleData& indelSampleData(indelDataPtr->getSampleData(sampleIndex));
                const auto& errorRates(indelSampleData.getErrorRates());
                const bool isIndelInCandAlignment(indelsInCandAlignment.count(evaluationIndel)!=0);
                if (isIndelInCandAlignment)
                {
                    // [[SAMPLE haplotype evaluationIndel TO OBSERVED haplotype (evaluationIndel+basecalling error)]]
                    //
                    // this represents the score of the read arising from the evaluationIndel haplotype, as explained by a gapless alignment:
                    //
                    updateIndelScoringInfo(indelScoringInfo,evaluationIndel,isIndelInCandAlignment,evaluationIndel,score,&candAlignment);

                    // [[SAMPLE haplotype reference TO OBSERVED haplotype (evaluationIndel+basecalling error)]]
                    //
                    // this represents the score of the read arising from the reference haplotype, as explained by a single indel error
                    //
                    updateIndelScoringInfo(indelScoringInfo,evaluationIndel,false,evaluationIndel,score+errorRates.refToIndelErrorProb.getLogValue(),&candAlignment);

                    // mark this as an alternate indel score for orthogonal indels:
                    for (const IndelKey& orthogonalIndel : orthogonalIndelMap[evaluationIndel])
                    {
                        // [[SAMPLE haplotype reference TO OBSERVED haplotype (evaluationIndel+basecalling error)]]
                        //
                        // same number as above, but provided in a format so that it is easily accessed from orthogonalIndel
                        //
                        updateIndelScoringInfo(indelScoringInfo, orthogonalIndel, false, orthogonalIndel,
                                               score + errorRates.refToIndelErrorProb.getLogValue(), &candAlignment);

                        {
                            // Note that this entry is designed to support the legacy "alt_indels" enumeration attached to ReadScoreData below, once that system is retired
                            //  (because indel overlap will be explicitely modeled), then this type of entry can be removed and the scoring structure can be simplified:

                            // [[SAMPLE haplotype evaluationIndel TO OBSERVED haplotype (evaluationIndel+basecalling error)]]
                            //
                            // same number as above, but provided in a format so that it is easily accessed from orthogonalIndel
                            //
                            updateIndelScoringInfo(indelScoringInfo,orthogonalIndel,true,evaluationIndel,score,&candAlignment);
                        }
                    }
                }
                else
                {
                    // if indel is not present, we must determine whether an orthogonal indel is present in the
                    // alignment to score this alignment in the right category:
                    //
                    const auto othogonalIndelIter(which_interfering_indel(indelsInCandAlignment,evaluationIndel));
                    const bool isIndelOrthogonalToCandAlignment(othogonalIndelIter != indelsInCandAlignment.end());

                    // record nonCandidate indels which are orthogonal to the evaluation set
                    // these will be handled separately below
                    //
                    // TODO what if more than one orthogonal indels are found?
                    if (isIndelOrthogonalToCandAlignment)
                    {
                        if (! indelsToEvaluate.count(*othogonalIndelIter))
                        {
                            nonCandidateIndelsOrthogonalToEvaluationIndels.insert(*othogonalIndelIter);
                        }
                    }

                    if (! isIndelOrthogonalToCandAlignment)
                    {
                        // [[SAMPLE haplotype reference TO OBSERVED haplotype (reference+basecalling error)]]
                        //
                        // this represents the score of the read arising from the reference haplotype
                        //
                        // note that "reference haplotype" means the reference allele over the reference span
                        // of the evaluation indel only
                        //
                        updateIndelScoringInfo(indelScoringInfo,evaluationIndel,false,evaluationIndel,score,&candAlignment);

                        // [[SAMPLE haplotype evaluationIndel TO OBSERVED haplotype (reference+basecalling error)]]
                        //
                        // this represents the score of the read arising from the evaluationIndel haplotype
                        // but with a single indel error reverting it to reference:
                        //
                        updateIndelScoringInfo(indelScoringInfo,evaluationIndel,true,evaluationIndel,score+errorRates.indelToRefErrorProb.getLogValue(),&candAlignment);
                    }
                    else
                    {
                        // [[SAMPLE haplotype evaluationIndel TO OBSERVED haplotype (orthogonalIndel+basecalling error)]]
                        //
                        // this represents the score of the read arising from the evaluationIndel haplotype
                        // but with a single indel error leading to an observation of the orthogonal haplotype
                        //
                        // note in the absence of and indel->indel2 error rate, we use the indel->ref rate as a proxy
                        //
                        updateIndelScoringInfo(indelScoringInfo,evaluationIndel,true,evaluationIndel,score+errorRates.indelToRefErrorProb.getLogValue(),&candAlignment);
                    }
                }
            }

            // run a separate loop to handle the remaining non-candidate indels in this alignment
            // which are orthogonal to at least one evaluation indel:
            //
            for (const IndelKey& nonCandidateIndel : nonCandidateIndelsOrthogonalToEvaluationIndels)
            {
                const IndelData* indelDataPtr(indelBuffer.getIndelDataPtr(nonCandidateIndel));
                assert(indelDataPtr != nullptr);
                const IndelSampleData& indelSampleData(indelDataPtr->getSampleData(sampleIndex));

                const auto& errorRates(indelSampleData.getErrorRates());

                // test whether indel is orthogonal to each evaluationIndel, if it is, then:
                //
                // record the (nonCandidateIndel from ref) score as a possible ref obs for the evaluationIndel
                //
                // ?TODO? -- also could make an observation on the (evalautionIndel TO nonCandidateIndel) score here,
                //          although it is handled above as well

                for (const IndelKey& evaluationIndel : indelsToEvaluate)
                {
                    const bool isIndelOrthogonalToCandAlignment(is_indel_conflict(nonCandidateIndel,evaluationIndel));
                    if (! isIndelOrthogonalToCandAlignment) continue;

                    // [[SAMPLE haplotype reference TO OBSERVED haplotype (nonCandidateIndel+basecalling error)]]
                    //
                    // this represents the score of the read arising from the reference haplotype
                    // but with a single indel error leading to an observation of the noncandidateIndel
                    //
                    updateIndelScoringInfo(indelScoringInfo,evaluationIndel,false,evaluationIndel,score+errorRates.refToIndelErrorProb.getLogValue(),&candAlignment);
                }
            }
        }
    }


    // (3) assemble indel calling scores for each indel into the indel buffer:
    //
    //
    const unsigned read_length(rseg.read_size());

    // for indel caller calculation we always need the "rest-of-genome"
    // alignment score, which is calculated from the number of non-ambiguous sites:
    uint16_t nonAmbiguousBasesInRead(0);
    {
        const bam_seq bseq(rseg.get_bam_read());
        for (unsigned i(0); i<read_length; ++i)
        {
            if (bseq.get_code(i) == BAM_BASE::ANY) continue;
            nonAmbiguousBasesInRead++;
        }
    }

    {
        const bool is_tier1_read(rseg.is_tier1_mapping());

        for (const IndelKey& evaluationIndel : indelsToEvaluate)
        {
            // we test for presence of the indel on the highest
            // scoring alignment because breakpoint overlap has
            // already been tested there, allowing us to skip this
            // step:
            const bool isIndelInMaxCandAlignment(indelsInMaxCandAlignment.count(evaluationIndel)!=0);

#ifdef DEBUG_ALIGN
            log_os << "VARMIT: final indel scan: " << evaluationIndel;
            log_os << "VARMIT: is_present_on_max_path?: " << isIndelInMaxCandAlignment << "\n";
#endif

            // sanity check that these exist in iks (unless search is incomplete):
            // 1) indel present
            // 2) indel absent
            // 3) conflicting indels present
            //

            // 1) indel present:
            double evaluationIndelScore(maxCandAlignmentScore);
            if (! isIndelInMaxCandAlignment)
            {
                const iks_map_t::iterator indelScoringIter(indelScoringInfo.find(std::make_pair(evaluationIndel,std::make_pair(true,evaluationIndel))));
                const bool is_found(indelScoringIter!=indelScoringInfo.end());
                if (! is_found)
                {
                    if (is_incomplete_search) continue;

                    if (is_safe_mode)
                    {
                        log_os << "WARNING: ";
                    }
                    else
                    {
                        log_os << "ERROR: ";
                    }

                    log_os << "failed to find expected alignment for indel: " << evaluationIndel
                           << "\twhile evaluating read_segment:\n" << rseg << "\n";

                    if (is_safe_mode)
                    {
                        continue;
                    }
                    else
                    {
                        exit(EXIT_FAILURE);
                    }
                }

                const iks_map_t::mapped_type& indelAlignmentInfo(indelScoringIter->second);
                evaluationIndelScore=indelAlignmentInfo.first;
            }

#ifdef DEBUG_ALIGN
            log_os << "VARMIT: called_indel_path_lnp " << evaluationIndelScore << "\n";
#endif

            // 2) indel absent w/o interference:
            double referenceScore(0);
            {
                const iks_map_t::iterator indelScoringIter(indelScoringInfo.find(std::make_pair(evaluationIndel,std::make_pair(false,evaluationIndel))));
                const bool is_found(indelScoringIter!=indelScoringInfo.end());
                if (is_incomplete_search && (! is_found)) continue;
                if (! is_found)
                {
                    if (is_safe_mode)
                    {
                        log_os << "WARNING: ";
                    }
                    else
                    {
                        log_os << "ERROR: ";
                    }

                    log_os << "failed to find reference alignment while evaluating read_segment:\n" << rseg << "\n";

                    if (is_safe_mode)
                    {
                        continue;
                    }
                    else
                    {
                        exit(EXIT_FAILURE);
                    }
                }

                referenceScore=indelScoringIter->second.first;

#ifdef DEBUG_ALIGN
                log_os << "VARMIT: ref_path_lnp " << referenceScore << "\n";
#endif
            }

            // assemble basic score data:
            //
            // where is the indel in fwd_strand read coordinates wrt the most likely read alignment?
            //
            /// TODO: should this be the most-likely alignment (present setting) or the most likely
            ///       alignment which contains the indel? If the former, we need to do something
            ///       about handling negative positions.
            ///
            //
            // To make this work, the reference range is expanded by one. This is because
            // getLowestFwdReadPosForRefRange needs to translate to  mapped read position,
            // not a position with a deletion. So it is effectively the mapped position adjacent
            // to the indel.
            const pos_t readPos(getLowestFwdReadPosForRefRange(maxCandAlignment.al,
                                                               known_pos_range(evaluationIndel.pos-1,
                                                                               evaluationIndel.right_pos()+1)) );

#ifdef DEBUG_ALIGN
            //            correct for strandedness
            log_os << "indelpos: " << evaluationIndel.pos << " readpos: + " << rseg.genome_align().pos << " rp: " << readpos << "\n";
#endif

            const unsigned fullReadLength(rseg.full_read_size());

            // this value is required for an RNA-Seq EVS feature
            pos_t distanceFromClosestReadEdge(fullReadLength);
            {
                alignment revMaxAlignment = maxCandAlignment.al;
                revMaxAlignment.is_fwd_strand = (not revMaxAlignment.is_fwd_strand);
                const pos_t revReadPos(getLowestFwdReadPosForRefRange(revMaxAlignment,
                                                                      known_pos_range(evaluationIndel.pos-1,
                                                                                      evaluationIndel.right_pos()+1)) );

                const unsigned fullReadOffset(rseg.full_read_offset());

                // care needs to be taken here, we want min(fullReadPos,fullRevReadPos), but if there is an error in
                // locating the indel's read pos from either direction, readPos or revReadPos may be set to -1.
                //
                if (readPos >= 0)
                {
                    distanceFromClosestReadEdge = (readPos + fullReadOffset);

                }
                if (revReadPos >= 0)
                {
                    const pos_t fullRevReadPos(revReadPos + (fullReadLength-(fullReadOffset+read_length)));
                    if (fullRevReadPos < distanceFromClosestReadEdge)
                    {
                        distanceFromClosestReadEdge = fullRevReadPos;
                    }
                }
            }

            ReadPathScores rps(referenceScore,evaluationIndelScore,nonAmbiguousBasesInRead,read_length,is_tier1_read,maxCandAlignment.al.is_fwd_strand,
                               (int16_t) readPos, (int16_t) distanceFromClosestReadEdge);

            // start adding alternate indel alleles, if present:

            for (const IndelKey& orthogonalIndel : orthogonalIndelMap[evaluationIndel])
            {
                const iks_map_t::iterator indelScoringIter(indelScoringInfo.find(std::make_pair(evaluationIndel,std::make_pair(true,orthogonalIndel))));

#ifdef DEBUG_ALIGN
                log_os << "VARMIT: alternate_indel " << orthogonalIndel;
#endif

                // TODO consider a way that the basic stores
                // could still be used even if we can't get all
                // alternate paths?
                //
                const bool is_found(indelScoringIter!=indelScoringInfo.end());

                //if(is_incomplete_search and (not is_found)) continue;
                //assert(is_found);

                // for now, this case has to be allowed, because it is possible
                // for no alignment not including the called indel to be in range
                // of the alternate indel
                //
                // for now we continue when the alt case is not found.
                //
                // TODO tighten the re-alignment procedure so that we know when the
                // alt case should be missing -- possibly use ref_path instead of continuing
                // in this case
                //
                if (is_found)
                {
                    rps.insertAlt(orthogonalIndel, indelScoringIter->second.first);
                    //                        rps.alt_indel[*k] = j->second.first;
                }
                else
                {
                    // TODO look at these cases in more detail before allowing them!
                    continue;
                    //rp.alt_indel[*k] = ref_path_lnp;
                }

#ifdef DEBUG_ALIGN
                //log_os << "VARMIT: alternate_indel_path_found: " << is_found << "\n";
                log_os << "VARMIT: alternate_indel_path_lnp " << indelScoringIter->second.first << "\n";
#endif
            }

            {
                IndelData* evaluationIndelDataPtr(indelBuffer.getIndelDataPtr(evaluationIndel));
                assert(nullptr != evaluationIndelDataPtr);
                IndelSampleData& evaluationIndelSampleData(evaluationIndelDataPtr->getSampleData(sampleIndex));
                evaluationIndelSampleData.read_path_lnp[rseg.getReadIndex()] = rps;
            }

#ifdef DEBUG_ALIGN
            log_os << "VARMIT: modified indel data: " << *(evaluationIndelDataPtr);
#endif
        }
    }
}


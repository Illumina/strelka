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

/// \file
/// \author Sangtae Kim
///

#include "ActiveRegionProcessor.hh"

#include "assembly/IterativeAssembler.hh"
#include "blt_util/algo_util.hh"

#include "boost/algorithm/string.hpp"
#include "boost/make_unique.hpp"

// compile with this macro to get verbose output:
//#define DEBUG_ACTIVE_REGION

void ActiveRegionProcessor::addHaplotypesToExclude(const std::vector<std::string>& haplotypeToExclude)
{
    _haplotypesToExclude = haplotypeToExclude;
}

const std::vector<std::string>& ActiveRegionProcessor::getSelectedHaplotypes() const
{
    return _selectedHaplotypes;
}

void ActiveRegionProcessor::processHaplotypes()
{
    // Check whether the active region is included in the read buffer
    const bool isRangeValid = (_posRange.begin_pos() >= _readBuffer.getBeginPos())
                              && (_posRange.end_pos() <= _readBuffer.getEndPos());

    // if the active region is not included in the read buffer or if it is too large,
    // bypass haplotyping
    if ((! isRangeValid) || (_posRange.size() > MaxRefSpanToBypassAssembly))
    {
        doNotUseHaplotyping();
    }
    else
    {
        bool isHaplotypingSuccess = generateHaplotypesWithCounting();
        if (not isHaplotypingSuccess)
        {
            // counting failed. Try assembly.
            isHaplotypingSuccess = generateHaplotypesWithAssembly();
        }

        if (isHaplotypingSuccess)
        {
            processSelectedHaplotypes();
        }
        else
        {
            // both counting and assembly failed
            // do not use haplotyping to determine indel candidacy
            doNotUseHaplotyping();
        }
    }
}

bool ActiveRegionProcessor::generateHaplotypesWithCounting()
{
    static const bool includePartialReads(false);
    ActiveRegionReadInfo readInfo;
    _readBuffer.getReadSegments(_posRange, readInfo, includePartialReads);

    // give up counting in the degenerate case of no read support:
    if (readInfo.numReadsAlignedToActiveRegion == 0) return false;

    const unsigned numReadsCoveringFullRegion((unsigned int) readInfo.readSegmentsForHaplotypeGeneration.size());

    // give up counting if there are not enough reads fully covering the region
    if (numReadsCoveringFullRegion < (MinFracReadsCoveringRegion*readInfo.numReadsAlignedToActiveRegion))
        return false;

    _numReadsUsedToGenerateHaplotypes = readInfo.numReadsAlignedToActiveRegion;
    HaplotypeToAlignIdSet haplotypeToAlignIdSet;
    for (const auto& entry : readInfo.readSegmentsForHaplotypeGeneration)
    {
        align_id_t alignId = entry.first;

        const std::string& haplotype(entry.second);

        if (!haplotypeToAlignIdSet.count(haplotype))
            haplotypeToAlignIdSet[haplotype] = std::vector<align_id_t>();
        haplotypeToAlignIdSet[haplotype].push_back(alignId);
    }

#ifdef DEBUG_ACTIVE_REGION
    std::cerr << _sampleIndex << "\t" << _posRange.begin_pos()+1 << '\t' << _posRange.end_pos() << '\t' << _refSegment << "\tCounting"<< std::endl;
#endif

    selectHaplotypes(haplotypeToAlignIdSet);

    return true;
}


bool ActiveRegionProcessor::generateHaplotypesWithAssembly()
{
    // Expand the region to include left/right anchors.
    // TODO: anchors may be too short if there are SNVs close to anchors
    // prefix anchor
    pos_t minBeginPos(0u);
    if (_posRange.begin_pos() > ActiveRegionReadBuffer::MaxAssemblyPadding)
        minBeginPos = _posRange.begin_pos() - ActiveRegionReadBuffer::MaxAssemblyPadding;
    minBeginPos = std::max(_readBuffer.getBeginPos(), minBeginPos);

    pos_t beginPos(_posRange.begin_pos());
    for (; beginPos > minBeginPos; --beginPos)
    {
        // anchor should not include a variant position
        if (_readBuffer.isCandidateVariant(beginPos-1)) break;
    }

    // suffix anchor
    pos_t maxEndPos = std::min(_readBuffer.getEndPos(), _posRange.end_pos() + ActiveRegionReadBuffer::MaxAssemblyPadding);
    pos_t endPos(_posRange.end_pos());
    for (; endPos < maxEndPos; ++endPos)
    {
        if (_readBuffer.isCandidateVariant(endPos)) break;
    }

    // prefix anchor ends with the first base of the active region
    std::string prefixAnchor;
    _ref.get_substring(beginPos, _posRange.begin_pos() - beginPos + 1, prefixAnchor);

    // suffix anchor starts with the last base of the active region
    std::string suffixAnchor;
    _ref.get_substring(_posRange.end_pos()-1, endPos-_posRange.end_pos() + 1, suffixAnchor);

    unsigned minReadSegmentLength((unsigned int) (prefixAnchor.size() + suffixAnchor.size()));

    // get read segments
    static const bool includePartialReads(true);
    ActiveRegionReadInfo readInfo;
    _readBuffer.getReadSegments(known_pos_range2(beginPos, endPos), readInfo, includePartialReads, minReadSegmentLength);

    /// \TODO: In the check below, it may be more consistent to replace "numReadsAlignedToActiveRegion" with the count
    ///        of reads which are eligible to go into the assembler: "readSegmentsForHaplotypeGeneration.size()"

    // too many reads; do not perform assembly (too time-consuming)
    if (readInfo.numReadsAlignedToActiveRegion > MinNumReadsToBypassAssembly)
        return false;   // assembly fail; bypass indels later

    AssemblyReadInput reads;
    std::vector<align_id_t> readIndexToAlignId;

    for (const auto& entry : readInfo.readSegmentsForHaplotypeGeneration)
    {
        align_id_t alignId = entry.first;

        const std::string& readSegment(entry.second);

        if (not readSegment.empty())
        {
            reads.push_back(readSegment);
            readIndexToAlignId.push_back(alignId);
        }
    }

    AssemblyReadOutput assemblyReadOutput;
    Assembly contigs;

    IterativeAssemblerOptions assembleOption;

    // We only accept haplotypes that start with prefixAnchor and end with suffixAnchor.
    // So, the minimum size of legitimate haplotypes is prefixAnchor.size() + suffixAnchor.size().
    // For most regions, minWordLength will be 20.
    // If there are SNVs closer to anchors, minWordLength may be smaller than 20.
    assembleOption.minWordLength = minReadSegmentLength;

    // maxWordLength must not be smaller than minWordLength.
    unsigned maxWordLength(std::max(minReadSegmentLength, ActiveRegionProcessor::MaxAssemblyWordSize));
    assembleOption.maxWordLength = maxWordLength;
    assembleOption.minCoverage = MinAssemblyCoverage;

    // perform assembly
    runIterativeAssembler(assembleOption, reads, assemblyReadOutput, contigs);

    unsigned totalNumReadsUsedInAssembly(0);
    for (const auto& assemblyReadInfo : assemblyReadOutput)
    {
        if (assemblyReadInfo.isUsed and (not assemblyReadInfo.isPseudo))
            ++totalNumReadsUsedInAssembly;
    }
    _numReadsUsedToGenerateHaplotypes = totalNumReadsUsedInAssembly;

    HaplotypeToAlignIdSet haplotypeToAlignIdSet;
    bool isNonRefHaplotypeFound(false);

#ifdef DEBUG_ACTIVE_REGION
    std::cerr << _sampleIndex << "\t" << _posRange.begin_pos()+1 << '\t' << _posRange.end_pos() << '\t' << _refSegment << "\tAssembly"<< std::endl;
#endif

    for (unsigned contigIndex(0); contigIndex<contigs.size(); ++contigIndex)
    {
        const std::string& contig(contigs[contigIndex].seq);
        // ignore if the contig does not contain prefix anchor
        auto start(contig.find(prefixAnchor));
        if (start == std::string::npos) continue;

        // remove prefix padding
        start += prefixAnchor.length() - 1;

        // ignore if the contig does not contain suffix anchor
        auto end(contig.rfind(suffixAnchor));
        if (end == std::string::npos or start > end) continue;

        // remove suffix padding
        end += 1;

        const std::string haplotype(contig.substr(start, end-start));

        auto alignIds = std::vector<align_id_t>();
        bool containsUniqueRead(false);
        for (unsigned readIndex : contigs[contigIndex].supportReads)
        {
            const auto& assemblyReadInfo(assemblyReadOutput[readIndex]);
            if (assemblyReadInfo.isPseudo)
            {
                // pseudo reads are ignored
                continue;
            }
            if ((not containsUniqueRead) and (assemblyReadInfo.contigIds.size() == 1))
                containsUniqueRead = true;
            alignIds.push_back(readIndexToAlignId[readIndex]);
        }

        // ignore if there's no read uniquely supporting the contig
        if (not containsUniqueRead) continue;

        if (haplotype != _refSegment)
            isNonRefHaplotypeFound = true;
        haplotypeToAlignIdSet[haplotype] = alignIds;
    }

    // assembly fails if no alt haplotype is found
    if (not isNonRefHaplotypeFound)
        return false;

    selectHaplotypes(haplotypeToAlignIdSet);

    return true;
}

void ActiveRegionProcessor::doNotUseHaplotyping()
{
#ifdef DEBUG_ACTIVE_REGION
    std::cerr << _sampleIndex << "\t" << _posRange.begin_pos()+1 << '\t' << _posRange.end_pos() << "\tBypass"<< std::endl;
#endif

    assert(_posRange.end_pos() > _posRange.begin_pos());

    auto it(_indelBuffer.positionIterator(_posRange.begin_pos()));
    const auto it_end(_indelBuffer.positionIterator(_posRange.end_pos()));

    for (; it!=it_end; ++it)
    {
        const IndelKey& indelKey(it->first);

        if (indelKey.is_breakpoint()) continue;

        // do not endorse mismatches if haplotyping fails
        if (indelKey.isMismatch()) continue;

        IndelData& indelData(getIndelData(it));
        indelData.activeRegionId = _posRange.begin_pos();

        // indicate that this indel was not discovered through haplotyping in this sample
        IndelSampleData& indelSampleData(indelData.getSampleData(_sampleIndex));
        indelSampleData.isHaplotypingBypassed = true;
    }
}



/// \return true if haplotypes are the same length and have exactly one mismatch
static
bool
doHaplotypesMeetPhasingErrorCondition1(
    const std::string& hap1,
    const std::string& hap2)
{
    // 2.
    if ((hap1.length() == hap2.length()) and (hap1 != hap2))
    {
        const auto retval = std::mismatch(hap1.begin(), hap1.end(), hap2.begin());
        const auto retval2 = std::mismatch(retval.first + 1, hap1.end(), retval.second + 1);
        if (retval2.first == hap1.end())
        {
            return true;
        }
    }
    return false;
}


/// \return true if haplotype2 appears to be a sequencer phasing error 'echo' of haplotype 1
///
/// test specifically for a very clean sequencer phasing error appearing as haplotype2 when the
/// sample contains a single true homozygous haplotype which we've detected in haplotype1
///
/// The error is found by testing for:
///
/// 1. the top two haplotypes are the same length and differ by only one basecall
/// 2. the second haplotype is observed exclusively on a single strand.
/// 3. that basecall difference is found at the (begin|end) of a homopolymer track, changing the base to match
///    the hpol base in the hap containing only (rev|fwd) support.
/// 4. the corresponding hpol is at least minPhaseErrorHpolSize in length (this is the hpol in hap1, not reference)
///
static
bool
isFilterSecondHaplotypeAsSequencerPhasingNoise(
    const ActiveRegionReadBuffer& readBuffer,
    const HaplotypeToAlignIdSet& haplotypeToAlignIdSet,
    const std::string& hap1,
    const std::string& hap2)
{
    static const int minPhaseErrorHpolSize(10);

    // test condition 1
    //
    // arrange this filter to come first so that we don't have to identify dups or read strands for most cases:
    if (not doHaplotypesMeetPhasingErrorCondition1(hap1,hap2)) return false;


    // test condition 2
    //

    // 2a: get haplotype supporting read lists:
    const auto hap1MapIter(haplotypeToAlignIdSet.find(hap1));
    assert(hap1MapIter != haplotypeToAlignIdSet.end());
    const auto& hap1AlignIdList(hap1MapIter->second);

    const auto hap2MapIter(haplotypeToAlignIdSet.find(hap2));
    assert(hap2MapIter != haplotypeToAlignIdSet.end());
    const auto& hap2AlignIdList(hap2MapIter->second);

    // 2b: identify duplicate reads:
    const std::set<align_id_t> dups(
        getDuplicatesInSortedInput(std::begin(hap1AlignIdList),std::end(hap1AlignIdList),
                                   std::begin(hap2AlignIdList),std::end(hap2AlignIdList)));

    // 2c: get total and unique counts for haplotype2:
    const unsigned hap2Count(hap2AlignIdList.size());
    const unsigned hap2UniqueCount(hap2Count-dups.size());

    // 2d: identify stranded counts:

    /// \return unique fwd-strand counts
    auto getHaplotypeNonDupFwdCount = [&](const std::vector<align_id_t>& alignIdList)
    {
        unsigned fwdCount(0);
        for (const auto alignId : alignIdList)
        {
            if (dups.count(alignId) > 0) continue;
            if (readBuffer.getAlignInfo(alignId).isForwardStrand) fwdCount++;
        }
        return fwdCount;
    };

    const unsigned hap2UniqueFwdCount(getHaplotypeNonDupFwdCount(hap2AlignIdList));
    assert(hap2UniqueFwdCount <= hap2UniqueCount);

    // 2e: test condition 2
    if ((hap2UniqueFwdCount > 0) and (hap2UniqueFwdCount < hap2UniqueCount)) return false;


    // test conditions 3 and 4
    //

    const auto retval = std::mismatch(hap1.begin(), hap1.end(), hap2.begin());
    if (hap2UniqueFwdCount == 0)
    {
        auto hap2iter(retval.second);
        const char hap2base(*hap2iter);
        for (; hap2iter != hap2.end(); hap2iter++)
        {
            if (*hap2iter != hap2base) break;
        }
        return ((hap2iter - retval.second) > minPhaseErrorHpolSize);
    }
    else
    {
        auto hap2iter(retval.second);
        const char hap2base(*hap2iter);
        while (true)
        {
            if (*hap2iter != hap2base) break;
            if (hap2iter == hap2.begin()) break;
            hap2iter--;
        }
        return ((retval.second - hap2iter) > minPhaseErrorHpolSize);
    }
}

void ActiveRegionProcessor::selectHaplotypes(const HaplotypeToAlignIdSet& haplotypeToAlignIdSet)
{
    using namespace std;
    vector<pair<unsigned,string>> haplotypeAndCounts;
    for (const auto& entry : haplotypeToAlignIdSet)
    {
        const string& haplotype(entry.first);
        const unsigned count(entry.second.size());

        if (count < MinHaplotypeCount) continue;

        haplotypeAndCounts.push_back(make_pair(count, haplotype));
    }

    if (haplotypeAndCounts.empty()) return;

    // sort in reverse order of counts
    sort(haplotypeAndCounts.begin(), haplotypeAndCounts.end(),
         [] (const pair<unsigned,string>& left, const pair<unsigned,string>& right)
    {
        return left.first > right.first;
    });

    const string& topHaplotype(haplotypeAndCounts[0].second);

    unsigned prevCount(std::numeric_limits<unsigned>::max());
    bool isReferenceSelected(false);

    // scan haplotypes and in reverse order of counts
    // and put up to (_ploidy+1) haplotypes into _selectedHaplotypes
    // haplotypes with the same count are stored in haplotypesWithSameCount
    // and added or dropped together
    // _selectedHaplotypes cannot store more than 2 alt haplotypes
    // e.g if _ploidy == 2
    // (15, ref), (12, hap1), (12, hap2) => _selectedHaplotypes = [ref, hap1, hap2]
    // (15, ref), (12, hap1), (12, hap2), (12, hap3) => _selectedHaplotypes = [ref]
    // (15, hap1), (12, hap2), (12, hap2) => _selectedHaplotypes = [hap1]
    vector<string> haplotypesWithSameCount;
    for (unsigned i(0); i<haplotypeAndCounts.size(); ++i)
    {
        const unsigned count(haplotypeAndCounts[i].first);
        const string& haplotype(haplotypeAndCounts[i].second);

        if (count < prevCount)
        {
            selectOrDropHaplotypesWithSameCount(
                haplotypesWithSameCount, haplotypeToAlignIdSet, isReferenceSelected);
        }

        if (_selectedHaplotypes.size() >= _ploidy) break;

        // apply sequencer phasing noise filter
        if (!isFilterSecondHaplotypeAsSequencerPhasingNoise(
                _readBuffer, haplotypeToAlignIdSet, topHaplotype, haplotype))
        {
            haplotypesWithSameCount.push_back(haplotype);
            if (haplotype == _refSegment) isReferenceSelected = true;
        }

        prevCount = count;
    }

    // add remaining haplotypes in haplotypesWithSameCount to _selectedHaplotypes
    if (!haplotypesWithSameCount.empty())
    {
        selectOrDropHaplotypesWithSameCount(
            haplotypesWithSameCount, haplotypeToAlignIdSet, isReferenceSelected);
    }
}

void ActiveRegionProcessor::selectOrDropHaplotypesWithSameCount(
    std::vector<std::string>& haplotypesWithSameCount,
    const HaplotypeToAlignIdSet& haplotypeToAlignIdSet,
    const bool isReferenceSelected)
{
    using namespace std;

    // add haplotypes in haplotypesWithSameCount to _selectedHaplotypes
    const unsigned numHapsWithTheSameCount(haplotypesWithSameCount.size());
    if (numHapsWithTheSameCount > 0)
    {
        // rename
        unsigned numHapsAfterAddingAll = (_selectedHaplotypes.size() + numHapsWithTheSameCount);
        // allow selecting up to _ploidy haplotypes
        // or _ploidy+1 haplotypes in case of tying count and reference is selected
        if ((numHapsAfterAddingAll <= _ploidy)
            || ((numHapsAfterAddingAll == _ploidy+1) && isReferenceSelected))
        {
            for (const string& haplotypeInBuffer : haplotypesWithSameCount)
            {
#ifdef DEBUG_ACTIVE_REGION
                std::cerr << haplotypeInBuffer << '\t' << haplotypeToAlignIdSet.at(haplotypeInBuffer).size() << std::endl;
#endif
                _selectedHaplotypes.push_back(haplotypeInBuffer);
                const std::vector<align_id_t>& alignIdList(haplotypeToAlignIdSet.at(haplotypeInBuffer));
                _selectedAlignIdLists.push_back(alignIdList);
            }
            haplotypesWithSameCount.clear();
        }
    }
}

void ActiveRegionProcessor::processSelectedHaplotypes()
{
    unsigned selectedHaplotypeIndex(0);
    std::vector<std::vector<IndelKey>> discoveredIndelsAndMismatches;
    std::vector<unsigned> selectedAltHaplotypeIndices;

    bool isAnyIndel(false);
    std::set<IndelKey> mismatchSet;

    for (const std::string& haplotype : _selectedHaplotypes)
    {
        if (haplotype != _refSegment)
        {
            int numIndels;
            std::vector<IndelKey> indelsAndMismatches;
            discoverIndelsAndMismatches(selectedHaplotypeIndex, indelsAndMismatches, numIndels);

            for (const auto& indelKey : indelsAndMismatches)
            {
                if (indelKey.isMismatch())
                {
                    mismatchSet.insert(indelKey);
                }
                else
                {
                    isAnyIndel = true;
                }
            }

            discoveredIndelsAndMismatches.push_back(indelsAndMismatches);
            selectedAltHaplotypeIndices.push_back(selectedHaplotypeIndex);
        }
        ++selectedHaplotypeIndex;
    }

    bool doNotAddMismatchesToIndelBuffer(false);
    const unsigned totalNumMismatches(mismatchSet.size());
    if (!isAnyIndel || (totalNumMismatches > MaxNumMismatchesToAddToIndelBuffer))
    {
        // for efficiency, if there is no indel or if there are too many mismatches,
        // don't add mismatches to the indel buffer
        doNotAddMismatchesToIndelBuffer = true;
    }

    HaplotypeId haplotypeId(0);
    for (unsigned selectedAltHaplotypeIndex : selectedAltHaplotypeIndices)
    {
        // alt haplotype gets non-zero haplotypeId
        ++haplotypeId;
        const std::vector<IndelKey>& indelsAndMismatches(discoveredIndelsAndMismatches[haplotypeId-1]);
        processDiscoveredIndelsAndMismatches(
            selectedAltHaplotypeIndex, haplotypeId, indelsAndMismatches, doNotAddMismatchesToIndelBuffer);
    }
}

void ActiveRegionProcessor::discoverIndelsAndMismatches(
    const unsigned selectedHaplotypeIndex,
    std::vector<IndelKey>& discoveredIndelsAndMismatches,
    int& numIndels)
{
    const std::string& haplotypeSeq(_selectedHaplotypes[selectedHaplotypeIndex]);
    assert (haplotypeSeq != _refSegment);

    pos_t referencePos;
    AlignmentResult<int> result;
    referencePos = _posRange.begin_pos();

    // Note that the aligner already left-shifts indels, but extra left-shifting is implemented below
    // because of reference edge artifacts, for example:
    //
    // There are cases that the active region was not triggered at the right position.
    // E.g. ref: GTCGAT, AR: TCGAT, Hap: T[ATAT]CGAT. In this case, T->TATAT should be left-shifted.
    //
    _aligner.align(haplotypeSeq.cbegin(),haplotypeSeq.cend(),_refSegment.cbegin(),_refSegment.cend(),result);

    const ALIGNPATH::path_t& alignPath = result.align.apath;

    pos_t haplotypePosOffset = 0;
    if (result.align.beginPos > 0)
    {
        assert(false && "Unexpected alignment segment");
    }

    numIndels = 0;
    unsigned numVariants(0);
    for (const auto& pathSegment : alignPath)
    {
        const unsigned segmentLength = pathSegment.length;

        switch (pathSegment.type)
        {
        case ALIGNPATH::SEQ_MATCH:
            referencePos += segmentLength;
            haplotypePosOffset += segmentLength;
            break;
        case ALIGNPATH::SEQ_MISMATCH:
        {
            for (unsigned i(0); i<segmentLength; ++i)
            {
                // two ARs may share one base and the same SNV could be find in both ARs
                // this if enables it to ignore the second discovery of the same SNV
                if (referencePos >= _prevActiveRegionEnd)
                {
                    // create a fake indelKey
                    // these are used in read realignment but not scored
                    auto insertSeq(haplotypeSeq.substr(haplotypePosOffset, 1));
                    auto mismatchIndelKey(IndelKey(referencePos, INDEL::MISMATCH, 1, insertSeq.c_str()));
                    discoveredIndelsAndMismatches.push_back(mismatchIndelKey);
                }

                ++referencePos;
                ++haplotypePosOffset;
            }
            ++numVariants;
            break;
        }
        case ALIGNPATH::INSERT:
        {
            if (segmentLength <= _maxIndelSize)
            {
                // left-align insertion
                // the insertion can be moved left by 1 base
                // if the last base of insertSeq equals to prevBase
                // E.g. GT -> GT(ATAT) vs GT -> G(TATA)T
                pos_t insertPos(referencePos);
                auto insertSeq(haplotypeSeq.substr((unsigned long) haplotypePosOffset, segmentLength));
                char prevBase = _ref.get_base(insertPos-1);
                while (insertSeq.back() == prevBase)
                {
                    // move insertion 1 base to left
                    insertSeq = prevBase + insertSeq;
                    insertSeq.pop_back();

                    --insertPos;
                    prevBase = _ref.get_base(insertPos-1);
                }

                // Ignore indel if (1) left-shifting causes it to overlap the previous AR
                // or (2) it extends past the right edge of the AR
                if ((prevBase != 'N') && (insertPos >= _prevActiveRegionEnd) && (insertPos < _posRange.end_pos()))
                {
                    discoveredIndelsAndMismatches.push_back(IndelKey(insertPos, INDEL::INDEL, 0, insertSeq.c_str()));
                    ++numIndels;
                    ++numVariants;
                }
            }
            haplotypePosOffset += segmentLength;
            break;
        }
        case ALIGNPATH::DELETE:
        {
            if (segmentLength <= _maxIndelSize)
            {
                // left-align deletion
                // the deletion can be moved left by 1 base
                // if the last base of the deleted sequence (lastDeletionBase) equals to prevBase
                // E.g. GT(ATAT) -> GT vs G(TATA)T -> GT
                pos_t deletePos(referencePos);
                char prevBase = _ref.get_base(deletePos-1);
                char lastDeletionBase = _ref.get_base(deletePos + segmentLength - 1);
                while (lastDeletionBase == prevBase)
                {
                    // move deletion 1 base to left
                    --deletePos;
                    lastDeletionBase = _ref.get_base(deletePos + segmentLength - 1);
                    prevBase = _ref.get_base(deletePos-1);
                }

                // Ignore indel if (1) left-shifting causes it to overlap the previous AR
                // or (2) it extends past the right edge of the AR
                if ((prevBase != 'N') && (deletePos >= _prevActiveRegionEnd) && (deletePos < _posRange.end_pos()))
                {
                    discoveredIndelsAndMismatches.push_back(IndelKey(deletePos, INDEL::INDEL, segmentLength));
                    ++numIndels;
                    ++numVariants;
                }
            }
            referencePos += segmentLength;
            break;
        }
        case ALIGNPATH::SOFT_CLIP:
        {
            referencePos += segmentLength;
            break;
        }
        default:
            assert(false && "Unexpected alignment segment");
        }
    }
}

void
ActiveRegionProcessor::processDiscoveredIndelsAndMismatches(
    const unsigned selectedHaplotypeIndex,
    const HaplotypeId haplotypeId,
    const std::vector<IndelKey>& discoveredIndelsAndMismatches,
    const bool doNotAddMismatchesToIndelBuffer)
{
    const std::vector<align_id_t>& alignIdList(_selectedAlignIdLists[selectedHaplotypeIndex]);

    // alignIdList.size()/_numReadsUsedToGenerateHaplotypes gives a meaning read count support ratio
    const float altHaplotypeCountRatio(alignIdList.size()/static_cast<float>(_numReadsUsedToGenerateHaplotypes));

    for (const auto& discoveredVariantKey : discoveredIndelsAndMismatches)
    {
        if (discoveredVariantKey.isMismatch())
        {
            _candidateSnvBuffer.addCandidateSnv(
                _sampleIndex, discoveredVariantKey.pos,
                discoveredVariantKey.insertSequence[0],
                haplotypeId, altHaplotypeCountRatio);
        }

        if (doNotAddMismatchesToIndelBuffer && discoveredVariantKey.isMismatch()) continue;

        for (const auto alignId : alignIdList)
        {
            const auto& alignInfo(_readBuffer.getAlignInfo(alignId));
            IndelObservationData indelObservationData;
            indelObservationData.iat = alignInfo.indelAlignType;
            indelObservationData.id = alignId;
            _indelBuffer.addIndelObservation(alignInfo.sampleIndex, {discoveredVariantKey, indelObservationData});
        }
        IndelData* indelDataPtr(_indelBuffer.getIndelDataPtr(discoveredVariantKey));
        assert((indelDataPtr != nullptr) && "Missing indelData");

        // Record active region ID in indel data
        indelDataPtr->activeRegionId = _posRange.begin_pos();

        // Update sample-specific indel details:
        IndelSampleData& indelSampleData(indelDataPtr->getSampleData(_sampleIndex));

        indelSampleData.haplotypeId += haplotypeId;

        assert (indelSampleData.haplotypeId <= 3);

        // All retios have the same denominator, so addition is valid:
        indelSampleData.altAlleleHaplotypeCountRatio += altHaplotypeCountRatio;
    }
}

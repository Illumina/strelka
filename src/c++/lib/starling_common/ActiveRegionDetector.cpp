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

#include "ActiveRegionDetector.hh"

ActiveRegionDetector::ActiveRegionDetector(
    const reference_contig_segment& ref,
    IndelBuffer& indelBuffer,
    CandidateSnvBuffer& candidateSnvBuffer,
    const unsigned maxIndelSize,
    const unsigned sampleCount,
    const bool isSomatic,
    const unsigned defaultPloidy) :
    _ref(ref),
    _sampleCount(sampleCount),
    _sampleActiveRegionDetector(sampleCount),
    _indelBuffer(indelBuffer),
    _candidateSnvBuffer(candidateSnvBuffer),
    _maxIndelSize(maxIndelSize),
    _isSomatic(isSomatic),
    _aligner(AlignmentScores<int>(ScoreMatch, ScoreMismatch, ScoreOpen, ScoreExtend, ScoreOffEdge, ScoreOpen, true, true))

{
    for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
    {
//        float minAlternativeAlleleFraction = MinAlternativeAlleleFraction;
//        if (_isSomatic && (sampleIndex != 0))
//            minAlternativeAlleleFraction = 0.1f;
        _sampleActiveRegionDetector[sampleIndex].reset(
            new SampleActiveRegionDetector(ref, MinAlternativeAlleleFraction, defaultPloidy, indelBuffer));
    }
}

ActiveRegionReadBuffer&
ActiveRegionDetector::getReadBuffer(const unsigned sampleIndex)
{
    assert (sampleIndex < _sampleCount);
    return _sampleActiveRegionDetector[sampleIndex]->getReadBuffer();
}

ActiveRegionId
ActiveRegionDetector::getActiveRegionId(const pos_t pos) const
{
    static const ActiveRegionId defaultActiveRegionId(-1);
    return _posToActiveRegionIdMap.getConstRefDefault(pos, defaultActiveRegionId);
}

void
ActiveRegionDetector::clearReadBuffer(const pos_t pos)
{
    for (unsigned sampleIndex(0); sampleIndex<_sampleCount; ++sampleIndex)
    {
        _sampleActiveRegionDetector[sampleIndex]->clearReadBuffer(pos);
    }
}

void
ActiveRegionDetector::updateEndPosition(const pos_t pos)
{
    processExistingActiveRegion(pos);
    for (unsigned sampleIndex(0); sampleIndex<_sampleCount; ++sampleIndex)
    {
        const auto activeRegionPtr = _sampleActiveRegionDetector[sampleIndex]->updateEndPosition(pos);
        if (activeRegionPtr)
        {
            updateActiveRegionRange(*activeRegionPtr);
        }
    }
}

void
ActiveRegionDetector::updateSamplePloidy(
    const unsigned sampleIndex, const pos_t pos, const unsigned ploidy)
{
    _sampleActiveRegionDetector[sampleIndex]->updatePloidy(pos, ploidy);
}

void
ActiveRegionDetector::clear()
{
    for (unsigned sampleIndex(0); sampleIndex<_sampleCount; ++sampleIndex)
    {
        const auto activeRegionPtr = _sampleActiveRegionDetector[sampleIndex]->closeActiveRegionDetector();
        if (activeRegionPtr)
        {
            updateActiveRegionRange(*activeRegionPtr);
        }
    }
    if (_synchronizedActiveRegion.end_pos())
        closeActiveRegion();
}

void
ActiveRegionDetector::clearPosToActiveRegionIdMapUpToPos(const pos_t pos)
{
    _posToActiveRegionIdMap.eraseTo(pos);
}

void
ActiveRegionDetector::updateActiveRegionRange(const ActiveRegion& sampleActiveRegion)
{
    if (_synchronizedActiveRegion.end_pos())
    {
        // _activeRegionRange is valid
        _synchronizedActiveRegion.merge_range(sampleActiveRegion);
    }
    else
    {
        _synchronizedActiveRegion = sampleActiveRegion;
    }
}

void
ActiveRegionDetector::processExistingActiveRegion(const pos_t pos)
{
    auto activeRegionEndPos(_synchronizedActiveRegion.end_pos());
    if (activeRegionEndPos &&
        (activeRegionEndPos <= (pos_t)(pos-ActiveRegionReadBuffer::MaxAssemblyPadding)))
    {
        bool canCloseActiveRegion(true);
        for (unsigned sampleIndex(0); sampleIndex<_sampleCount; ++sampleIndex)
        {
            const auto unclosedActiveRegionStartPos(getSampleActiveRegionDetector(sampleIndex).getActiveRegionStartPos());
            if ((unclosedActiveRegionStartPos >= 0) && (unclosedActiveRegionStartPos < activeRegionEndPos-1))
            {
                canCloseActiveRegion = false;
                break;
            }
        }

        if (canCloseActiveRegion)
            closeActiveRegion();
    }
}

void
ActiveRegionDetector::filterOutConflictingForcedIndels()
{
    // first check if one or more candidate indel is discovered
    // in this active region
    bool isCandidateIndel(false);

    auto it(_indelBuffer.positionIterator(_synchronizedActiveRegion.begin_pos()));
    auto it_end(_indelBuffer.positionIterator(_synchronizedActiveRegion.end_pos()));
    for (; it!=it_end; ++it)
    {
        IndelData& indelData(getIndelData(it));
        if (!indelData.isDiscoveredInActiveRegion()) continue;

        const IndelKey& indelKey(it->first);

        if (indelKey.isMismatch()) continue;

        bool isThisIndelCandidate = _indelBuffer.isCandidateIndel(indelKey, indelData);
        if (isThisIndelCandidate) isCandidateIndel = true;
    }

    if (! isCandidateIndel)
    {
        // if there exist no discovered candidate indel
        // there's no need to worry about forced indel conflicts
        return;
    }

    // find all the forced indels and set doNotGenotype=true
    it = _indelBuffer.positionIterator(_synchronizedActiveRegion.begin_pos());
    it_end = _indelBuffer.positionIterator(_synchronizedActiveRegion.end_pos());

    for (; it!=it_end; ++it)
    {
        IndelData& indelData(getIndelData(it));
        if (! indelData.isForcedOutput) continue;

        if (! indelData.isDiscoveredInActiveRegion())
        {
            // this forced indel was not discovered in AR
            // and may conflict with internal indels
            indelData.doNotGenotype = true;
        }
    }
}

void
ActiveRegionDetector::closeActiveRegion()
{
    assert (_synchronizedActiveRegion.size() > 0);

    std::vector<std::string> normalHaplotypes;

    for (unsigned sampleIndex(0); sampleIndex<_sampleCount; ++sampleIndex)
    {
        const unsigned ploidy(getPloidy(sampleIndex, _synchronizedActiveRegion));
        if (!_isSomatic)
        {
            ActiveRegionProcessor activeRegionProcessor(_synchronizedActiveRegion, _prevActiveRegionEnd,
                                                        _ref, _maxIndelSize, sampleIndex, ploidy,
                                                        _aligner, _sampleActiveRegionDetector[sampleIndex]->getReadBuffer(),
                                                        _indelBuffer, _candidateSnvBuffer);
            activeRegionProcessor.processHaplotypes();
        }
        else
        {
            // These are experimental and currently not enabled.
            if (sampleIndex == 0)
            {
                ActiveRegionProcessor activeRegionProcessor(_synchronizedActiveRegion, _prevActiveRegionEnd,
                                                            _ref, _maxIndelSize, sampleIndex, ploidy,
                                                            _aligner, _sampleActiveRegionDetector[sampleIndex]->getReadBuffer(),
                                                            _indelBuffer, _candidateSnvBuffer);
                activeRegionProcessor.processHaplotypes();
                normalHaplotypes = activeRegionProcessor.getSelectedHaplotypes();
            }
            else
            {
                ActiveRegionProcessor activeRegionProcessor(_synchronizedActiveRegion, _prevActiveRegionEnd,
                                                            _ref, _maxIndelSize, sampleIndex, 1u,
                                                            _aligner, _sampleActiveRegionDetector[sampleIndex]->getReadBuffer(),
                                                            _indelBuffer, _candidateSnvBuffer);
                activeRegionProcessor.addHaplotypesToExclude(normalHaplotypes);
                activeRegionProcessor.processHaplotypes();
            }
        }
    }

    filterOutConflictingForcedIndels();

    setPosToActiveRegionIdMap(_synchronizedActiveRegion);
    _prevActiveRegionEnd = _synchronizedActiveRegion.end_pos();
    _synchronizedActiveRegion.clear();
}

SampleActiveRegionDetector&
ActiveRegionDetector::getSampleActiveRegionDetector(unsigned sampleIndex)
{
    assert(sampleIndex < _sampleCount);
    return *(_sampleActiveRegionDetector[sampleIndex]);
}


unsigned
ActiveRegionDetector::getPloidy(const unsigned sampleIndex, const ActiveRegion activeRegion) const
{
    return _sampleActiveRegionDetector[sampleIndex]->getPloidy(activeRegion);
}


void
ActiveRegionDetector::setPosToActiveRegionIdMap(const ActiveRegion& activeRegionRange)
{
    if (activeRegionRange.size() > ActiveRegionProcessor::MaxRefSpanToBypassAssembly)
        return;

    ActiveRegionId activeRegionId(activeRegionRange.begin_pos());
    for (pos_t pos(activeRegionRange.begin_pos()); pos<activeRegionRange.end_pos(); ++pos)
    {
        _posToActiveRegionIdMap.getRef(pos) = activeRegionId;
    }
}

unsigned
SampleActiveRegionDetector::getPloidy(const pos_t pos)
{
    return _posToPloidyMap.getConstRefDefault(pos, _defaultPloidy);
}

unsigned
SampleActiveRegionDetector::getPloidy(const ActiveRegion& activeRegion)
{
    const auto leftPloidy(getPloidy(activeRegion.begin_pos()));
    const auto rightPloidy(getPloidy(activeRegion.end_pos()-1));

    return std::max(leftPloidy, rightPloidy);
}

void
SampleActiveRegionDetector::clearCoordinates()
{
    _prevAnchorPos = -1;
    _activeRegionStartPos = -1;
    _anchorPosFollowingPrevVariant = -1;
    _prevVariantPos = -1;
    _numVariants = 0;
}

void
SampleActiveRegionDetector::clearReadBuffer(const pos_t pos)
{
    _readBuffer.clearPos(pos);
    _posToPloidyMap.eraseTo(pos);
}



std::unique_ptr<ActiveRegion>
SampleActiveRegionDetector::createActiveRegion()
{
    // close the existing active region
    assert (_activeRegionStartPos < _anchorPosFollowingPrevVariant);

    auto activeRegionStartPos(_activeRegionStartPos);
    auto activeRegionEndPos(_anchorPosFollowingPrevVariant);

    // we have no existing active region at this point
    _numVariants = 0;
    _activeRegionStartPos = 0;

    return std::unique_ptr<ActiveRegion>(new ActiveRegion(activeRegionStartPos, activeRegionEndPos+1));
}

void
SampleActiveRegionDetector::updatePloidy(const pos_t pos, const unsigned ploidy)
{
    _posToPloidyMap.getRef(pos) = ploidy;
}

std::unique_ptr<ActiveRegion>
SampleActiveRegionDetector::updateEndPosition(const pos_t pos)
{
    _readBuffer.setEndPos(pos+1);

    if (_isBeginning)
    {
        _activeRegionStartPos = pos;
        _anchorPosFollowingPrevVariant = pos;
        _prevAnchorPos = pos;
        _isBeginning = false;
    }

    pos_t posToProcess(pos-1);
    if (posToProcess < 0) return nullptr;


    bool isCurrentPosCandidateVariant = _readBuffer.isCandidateVariant(posToProcess);
    const bool isDepthZero = _readBuffer.isDepthZero(posToProcess);
    // depth 0 position can be a candidate variant pos
    // but cannot open a new active region
    if (isDepthZero and (_numVariants == 0u))
    {
        isCurrentPosCandidateVariant = false;
    }

    const bool isAnchor = _readBuffer.isAnchor(posToProcess) and (not isCurrentPosCandidateVariant);

    if (!isCurrentPosCandidateVariant and !isAnchor) return nullptr;

//    std::cout << (posToProcess+1) << "\t" << isCurrentPosCandidateVariant << '\t' << isAnchor << std::endl;

    unsigned distanceFromPrevVariant = (unsigned) (posToProcess - _prevVariantPos);
    std::unique_ptr<ActiveRegion> activeRegion;
    if (distanceFromPrevVariant > MaxDistanceBetweenTwoVariants and _anchorPosFollowingPrevVariant >= 0)
    {
        if (_numVariants >= MinNumVariantsPerRegion )
        {
            activeRegion = createActiveRegion();
        }
        else
        {
            // no active region is detected
            _numVariants = 0;
        }
    }

    if (isAnchor)
    {
        if (_numVariants == 0)
        {
            // this position can be a start of an active region
            _activeRegionStartPos = posToProcess;
        }

        if (_anchorPosFollowingPrevVariant < 0)
        {
            _anchorPosFollowingPrevVariant = posToProcess;
        }
        _prevAnchorPos = posToProcess;
    }

    if (isCurrentPosCandidateVariant)
    {
        if (not _activeRegionStartPos)
            _activeRegionStartPos = _prevAnchorPos;
        // extend the existing active region
        ++_numVariants;
        _prevVariantPos = posToProcess;
        _anchorPosFollowingPrevVariant = -1;
    }

    return activeRegion;
}



std::unique_ptr<ActiveRegion>
SampleActiveRegionDetector::closeActiveRegionDetector()
{
    if (_isBeginning) return nullptr;

    std::unique_ptr<ActiveRegion> activeRegion;
    if (_numVariants >= MinNumVariantsPerRegion)
    {
        bool isStartPosUndetermined = (_activeRegionStartPos < 0);
        bool isEndPosUndetermined = (_anchorPosFollowingPrevVariant < 0);

        // if both the start and end pos are undetermined
        // and read buffer size is zero, do not close the active region
        if (!isStartPosUndetermined || !isEndPosUndetermined
            || (_readBuffer.getBeginPos() < _readBuffer.getEndPos()))
        {
            if (isStartPosUndetermined)
                _activeRegionStartPos = _readBuffer.getBeginPos();
            if (isEndPosUndetermined)
                _anchorPosFollowingPrevVariant = (_readBuffer.getEndPos()-1);
            activeRegion = createActiveRegion();
        }
    }

    clearCoordinates();

    return activeRegion;
}


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
/// \author Mitch Bekritsky
///

#include "IndelBuffer.hh"

#include "calibration/IndelErrorModel.hh"
#include "starling_common/AlleleReportInfoUtil.hh"

#include <iostream>



unsigned
IndelBuffer::
registerSample(
    const depth_buffer& db,
    const depth_buffer& db2,
    const bool isCountTowardsDepthFilter)
{
    assert(! _isFinalized);
    const unsigned sampleIndex(_indelSampleData.size());
    _indelSampleData.emplace_back(db,db2,isCountTowardsDepthFilter);
    return sampleIndex;
}



std::pair<IndelBuffer::iterator,IndelBuffer::iterator>
IndelBuffer::
rangeIterator(
    const pos_t begin_pos,
    const pos_t end_pos)
{
    const IndelKey end_range_key(end_pos);
    const iterator end(_indelBuffer.lower_bound(end_range_key));
    const IndelKey begin_range_key(begin_pos-static_cast<pos_t>(_opt.maxIndelSize));
    iterator begin(_indelBuffer.lower_bound(begin_range_key));
    for (; begin!=end; ++begin)
    {
        if (begin->first.right_pos() >= begin_pos) break;
    }
    return std::make_pair(begin,end);
}



// The goal is to return all indels with a left or right breakpoint in the
// range [begin_pos,end_pos]. Returning indels in addition to this set is
// acceptable.
//
// The indels_keys: "end_range_key" and "begin_range_key" take
// advantage of the indel NONE type, which sorts ahead of all other
// types at the same position.
//
std::pair<IndelBuffer::const_iterator,IndelBuffer::const_iterator>
IndelBuffer::
rangeIterator(
    const pos_t begin_pos,
    const pos_t end_pos) const
{
    const IndelKey end_range_key(end_pos);
    const const_iterator end(_indelBuffer.lower_bound(end_range_key));
    const IndelKey begin_range_key(begin_pos-static_cast<pos_t>(_opt.maxIndelSize));
    const_iterator begin(_indelBuffer.lower_bound(begin_range_key));
    for (; begin!=end; ++begin)
    {
        if (begin->first.right_pos() >= begin_pos) break;
    }
    return std::make_pair(begin,end);
}



bool
IndelBuffer::
addIndelObservation(
    const unsigned sampleIndex,
    const IndelObservation& obs)
{
    assert(_isFinalized);
    assert(obs.key.type != INDEL::NONE);

    // if not previously observed
    iterator indelIter(_indelBuffer.find(obs.key));
    const bool isNovel(indelIter == _indelBuffer.end());
    if (isNovel)
    {
        const auto retval = _indelBuffer.insert(std::make_pair(obs.key,IndelData(getSampleCount(), obs.key)));
        indelIter = retval.first;
    }

    IndelData& indelData(getIndelData(indelIter));
    if (isNovel)
    {
        indelData.initializeAuxInfo(_opt,_dopt, _ref);

        const auto& indelKey(indelIter->first);
        const bool isPrimitive = (
                                     indelKey.isMismatch() ||
                                     indelKey.isPrimitiveInsertionAllele() ||
                                     indelKey.isPrimitiveDeletionAllele());

        if (! isPrimitive)
        {
            // complex alleles are not genotyped
            indelData.doNotGenotype = true;
        }
    }
    indelData.addIndelObservation(sampleIndex, obs.data);

    return isNovel;
}



bool
IndelBuffer::
isCandidateIndelImplTestSignalNoise(
    const IndelKey& indelKey,
    const IndelData& indelData) const
{
    // determine if the signal/noise ratio of the indel qualifies it for candidacy in at least one sample
    //
    const unsigned sampleCount(indelData.getSampleCount());
    for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
    {
        const IndelSampleData& indelSampleData(indelData.getSampleData(sampleIndex));

        // for each sample, get the number of tier 1 reads supporting the indel
        // and the total number of tier 1 reads at this locus
        const unsigned tier1ReadSupportCount(indelSampleData.tier1_map_read_ids.size());
        unsigned totalReadCount(ebuff(sampleIndex).val(indelKey.pos-1));

        // total reads and indel reads are measured in different ways here, so the nonsensical
        // result of indel_reads>total_reads is possible. The total is fudged below to appear sane
        // before going into the count test:
        totalReadCount = std::max(totalReadCount,tier1ReadSupportCount);

        if (totalReadCount == 0) continue;

        // this min indel support coverage is applied regardless of error model:
        static const unsigned min_candidate_cov_floor(2);
        if (totalReadCount < min_candidate_cov_floor) continue;

#ifdef USE_SIMPLE_FREQ_FOR_INDEL_CANDIDACY
        // Test if the observed fraction of indel allele support is above a constant
        // threshold to determine indel candidacy
        //
        // This simplified backup candidacy model is maintained and often used for
        // study/optimization of indel error rates, so that impact of rates on quality
        // scores/genotype assignments and impact of rates on candidacy are not
        // confounded in a way which might otherwise make optimization concusions
        // difficult.
        //
        const double supportFraction(tier1ReadSupportCount/static_cast<double>(totalReadCount));
        if (supportFraction >= 0.10) return true;
#else
        // Test to see if the observed indel coverage has a binomial exact test
        // p-value above the rejection threshold. If this does not occur for the
        // counts observed in any sample, the indel cannot become a candidate
        if (_countCache.isRejectNull(totalReadCount, indelSampleData.getErrorRates().candidateIndelToRefErrorProb.getValue(), tier1ReadSupportCount))
        {
            return true;
        }
#endif
    }
    return false;
}



bool
IndelBuffer::
isCandidateIndelImplTestWeakSignal(
    const IndelData& indelData) const
{
    const unsigned sampleCount(indelData.getSampleCount());
    for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
    {
        const IndelSampleData& indelSampleData(indelData.getSampleData(sampleIndex));
        const unsigned tier1ReadSupportCount(indelSampleData.tier1_map_read_ids.size());
        static const unsigned min_candidate_cov_floor(1);
        if (tier1ReadSupportCount < min_candidate_cov_floor) continue;
        return true;
    }
    return false;
}



bool
IndelBuffer::
isCandidateIndelImplTest(
    const IndelKey& indelKey,
    const IndelData& indelData) const
{
    if (indelData.doNotGenotype) return false;

    // if haplotyping is enabled, indels not confirmed in active region are not candidate
    if (_opt.isHaplotypingEnabled && (!indelData.isDiscoveredInActiveRegion())) return false;

    if (_opt.is_candidate_indel_signal_test)
    {
        if (not isCandidateIndelImplTestSignalNoise(indelKey, indelData)) return false;
    }
    else
    {
        if (not isCandidateIndelImplTestWeakSignal(indelData)) return false;
    }

    //
    // test against short open-ended segments:
    //
    // call getInsertSize() instead of asking for the insert seq
    // so as to not finalize any incomplete insertions:
    if (indelKey.is_breakpoint() &&
        (_opt.min_candidate_indel_open_length > indelData.getBreakpointInsertSize()))
    {
        return false;
    }

    //
    // test against max_depth:
    //
    if (_maxCandidateDepth > 0.)
    {
        const pos_t depthPos(indelKey.pos-1);
        const unsigned sampleCount(getSampleCount());
        double estimatedLocusDepth(0);
        for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
        {
            if (not getIndelSampleData(sampleIndex).isCountTowardsDepthFilter) continue;

            const unsigned estdepth(ebuff(sampleIndex).val(depthPos));
            const unsigned estdepth2(ebuff2(sampleIndex).val(depthPos));
            estimatedLocusDepth += (estdepth + estdepth2);
        }
        if (estimatedLocusDepth > _maxCandidateDepth) return false;
    }

    return true;
}



void
IndelBuffer::
isCandidateIndelImpl(
    const IndelKey& indelKey,
    const IndelData& indelData) const
{
    bool isCandidate(isCandidateIndelImplTest(indelKey, indelData));

    // check whether the candidate has been externally specified:
    if (! isCandidate)
    {
        if ((! indelData.doNotGenotype) && indelData.is_external_candidate)
        {
            assert(indelKey.isPrimitiveDeletionAllele() || indelKey.isPrimitiveInsertionAllele());

            // primitive insertions and deletions are promoted to candidate status
            // even if there is no read support
            isCandidate = true;

            // mark that this indel doesn't have enough read support
            // even though it has the candidate status
            indelData.status.notDiscoveredFromReads = true;
        }
    }

    indelData.status.is_candidate_indel = isCandidate;
    indelData.status.is_candidate_indel_cached = true;
}



void
IndelBuffer::
findDataException(const IndelKey& indelKey) const
{
    std::ostringstream oss;
    oss << "Could not find indel_data for indel: " << indelKey;
    throw blt_exception(oss.str().c_str());
}



void
IndelBuffer::
clearIndelsAtPosition(const pos_t pos)
{
    const iterator i_begin(positionIterator(pos));
    const iterator i_end(positionIterator(pos + 1));
    _indelBuffer.erase(i_begin,i_end);
}



static
void
dump_range(
    IndelBuffer::const_iterator i,
    const IndelBuffer::const_iterator i_end,
    std::ostream& os)
{
    for (; i!=i_end; ++i)
    {
        os << "INDEL_KEY: " << i->first;
        os << "INDEL_DATA:\n";
        os << getIndelData(i);
    }
}



void
IndelBuffer::
dumpPosition(
    const pos_t pos,
    std::ostream& os) const
{
    dump_range(positionIterator(pos), positionIterator(pos + 1),os);
}



void
IndelBuffer::
dump(std::ostream& os) const
{
    os << "INDEL_BUFFER DUMP ON\n";
    dump_range(_indelBuffer.begin(),_indelBuffer.end(),os);
    os << "INDEL_BUFFER DUMP OFF\n";
}

void
IndelBuffer::
setMaxCandidateDepth(
    const double maxCandidateDepth)
{
    _maxCandidateDepth = maxCandidateDepth;
}

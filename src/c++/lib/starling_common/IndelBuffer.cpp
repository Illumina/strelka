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
/// \author Mitch Bekritsky
///

#include "IndelBuffer.hh"

#include "blt_util/prob_util.hh"
#include "calibration/IndelErrorModel.hh"

#include <iostream>



unsigned
IndelBuffer::
registerSample(
    const depth_buffer& db,
    const depth_buffer& db2,
    const double max_depth)
{
    assert(! _isFinalized);
    const unsigned sampleIndex(_indelSampleData.size());
    _indelSampleData.emplace_back(db,db2,max_depth);
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
    const IndelKey begin_range_key(begin_pos-static_cast<pos_t>(_opt.max_indel_size));
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
    const IndelKey begin_range_key(begin_pos-static_cast<pos_t>(_opt.max_indel_size));
    const_iterator begin(_indelBuffer.lower_bound(begin_range_key));
    for (; begin!=end; ++begin)
    {
        if (begin->first.right_pos() >= begin_pos) break;
    }
    return std::make_pair(begin,end);
}



static
void
scaleIndelErrorRate(
    const double logScaleFactor,
    double& indelErrorRate)
{
    static const double minIndelErrorProb(0.0);
    static const double maxIndelErrorProb(0.5);

    indelErrorRate = std::min(indelErrorRate, maxIndelErrorProb);
    indelErrorRate = softMaxInverseTransform(indelErrorRate, minIndelErrorProb, maxIndelErrorProb);
    indelErrorRate += logScaleFactor;
    indelErrorRate = softMaxTransform(indelErrorRate, minIndelErrorProb, maxIndelErrorProb);
}



bool
IndelBuffer::
addIndelObservation(
    const unsigned sampleId,
    const IndelObservation& obs)
{
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
    indelData.addIndelObservation(sampleId, obs.data);
    return isNovel;
}



bool
IndelBuffer::
isCandidateIndelImplTestSignalNoise(
    const IndelKey& indelKey,
    const IndelData& indelData) const
{
    // determine if the observed counts are sig wrt the error rate for at least one sample:
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

        // test to see if the observed indel coverage has a binomial exact test
        // p-value above the rejection threshold. If this does not occur for the
        // counts observed in any sample, the indel cannot become a candidate
        if (_countCache.isRejectNull(totalReadCount, indelData.errorRates.indelToRefErrorProb.getValue(), tier1ReadSupportCount))
        {
            return true;
        }
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
    // initialize indel error rates:
    if (! indelData.errorRates.isInit)
    {
        // get standard rates:
        starling_indel_report_info indelReportInfo;
        get_starling_indel_report_info(indelKey, indelData, _ref, indelReportInfo);
        double refToIndelErrorProb;
        double indelToRefErrorProb;
        _dopt.getIndelErrorModel().getIndelErrorRate(indelKey, indelReportInfo, refToIndelErrorProb, indelToRefErrorProb);
        indelData.errorRates.refToIndelErrorProb.updateValue(refToIndelErrorProb);
        indelData.errorRates.indelToRefErrorProb.updateValue(indelToRefErrorProb);

        // compute scaled rates:
        double scaledRefToIndelErrorProb = indelData.errorRates.refToIndelErrorProb.getValue();
        double scaledIndelToRefErrorProb = indelData.errorRates.indelToRefErrorProb.getValue();
        if (_opt.isIndelErrorRateFactor)
        {
            scaleIndelErrorRate(_dopt.logIndelErrorRateFactor, scaledRefToIndelErrorProb);
            scaleIndelErrorRate(_dopt.logIndelErrorRateFactor, scaledIndelToRefErrorProb);
        }

        indelData.errorRates.scaledRefToIndelErrorProb.updateValue(scaledRefToIndelErrorProb);
        indelData.errorRates.scaledIndelToRefErrorProb.updateValue(scaledIndelToRefErrorProb);

        indelData.errorRates.isInit = true;
    }

    // check whether the candidate has been externally specified:
    if (indelData.is_external_candidate) return true;

    if (_opt.is_candidate_indel_signal_test)
    {
        if (!isCandidateIndelImplTestSignalNoise(indelKey, indelData)) return false;
    }
    else
    {
        if (!isCandidateIndelImplTestWeakSignal(indelData)) return false;
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
    const unsigned sampleCount(getSampleCount());
    for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
    {
        const double max_depth(getIndelSampleData(sampleIndex).max_depth);
        if (max_depth <= 0.) continue;

        const unsigned estdepth(ebuff(sampleIndex).val(indelKey.pos-1));
        const unsigned estdepth2(ebuff2(sampleIndex).val(indelKey.pos-1));
        if ((estdepth+estdepth2) > max_depth) return false;
    }

    return true;
}



void
IndelBuffer::
isCandidateIndelImpl(
    const IndelKey& indelKey,
    const IndelData& indelData,
    const bool setIsCandidateIndel) const
{
    const bool is_candidate(isCandidateIndelImplTest(indelKey, indelData));
    if (setIsCandidateIndel)
        indelData.status.is_candidate_indel = is_candidate;
    indelData.status.is_candidate_indel_cached = true;
}



void
IndelBuffer::
findDataException(const IndelKey& indelKey) const
{
    std::ostringstream oss;
    oss << "ERROR: could not find indel_data for indel: " << indelKey;
    throw blt_exception(oss.str().c_str());
}



void
IndelBuffer::
clearPosition(const pos_t pos)
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

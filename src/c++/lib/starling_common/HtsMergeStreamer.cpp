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
/// \author Chris Saunders
///


#include "HtsMergeStreamer.hh"
#include "common/Exceptions.hh"

#include "boost/optional.hpp"



HtsMergeStreamer::
HtsMergeStreamer(
    const std::string& referenceFilename)
    : _referenceFilename(referenceFilename)
{}



const HtsMergeStreamer::OrderData&
HtsMergeStreamer::
getOrderForType(
    const unsigned orderIndex,
    const HTS_TYPE::index_t expectedHtsType) const
{
    const OrderData& orderData(_order[orderIndex]);
    assert(orderData.htsType == expectedHtsType);
    return orderData;
}



void
HtsMergeStreamer::
queueItem(
    const unsigned orderIndex)
{
    const auto htsType(getHtsType(orderIndex));

    boost::optional<pos_t> nextItemPos;
    if     (HTS_TYPE::BAM == htsType)
    {
        bam_streamer& bs(getHtsStreamer(orderIndex, _data._bam));
        if (bs.next())
        {
            nextItemPos = (bs.get_record_ptr()->pos() - 1);
        }
    }
    else if (HTS_TYPE::VCF == htsType)
    {
        vcf_streamer& vs(getHtsStreamer(orderIndex, _data._vcf));
        if (vs.next())
        {
            nextItemPos = (vs.get_record_ptr()->pos - 1);
        }
    }
    else if (HTS_TYPE::BED == htsType)
    {
        bed_streamer& bes(getHtsStreamer(orderIndex, _data._bed));
        if (bes.next())
        {
            nextItemPos = (bes.get_record_ptr()->begin - 1);
        }
    }

    if (nextItemPos)
    {
        _streamQueue.emplace(*nextItemPos, orderIndex);
    }
}



void
HtsMergeStreamer::
resetRegion(const std::string& region)
{
    assert(! region.empty());

    _region = region;
    _isStreamBegin = false;
    _isStreamEnd = false;
    _streamQueue = queue_t(); // why no .clear() for queues?

    const unsigned streamCount(_order.size());
    for (unsigned streamIndex(0); streamIndex < streamCount; ++streamIndex)
    {
        const auto& orderData(_order[streamIndex]);
        if (orderData.htsType == HTS_TYPE::BAM)
        {
            getHtsStreamer(streamIndex, _data._bam).resetRegion(region.c_str());
        }
        else if (orderData.htsType == HTS_TYPE::BED)
        {
            getHtsStreamer(streamIndex, _data._bed).resetRegion(region.c_str());
        }
        else if (orderData.htsType == HTS_TYPE::VCF)
        {
            getHtsStreamer(streamIndex, _data._vcf).resetRegion(region.c_str());
        }
        else
        {
            assert(false and "Unexpected hts file type.");
        }

        queueItem(streamIndex);
    }
}



bool
HtsMergeStreamer::
next()
{
    if (_isStreamEnd) return false;

    if (_isStreamBegin)
    {
        // reload stream_queue with current type and sample_no;
        queueNextItem();
        const HtsRecordSortData last = getCurrent();
        _streamQueue.pop();

        if (getCurrentPos() < last.pos)
        {
            using namespace illumina::common;
            const auto htsType(getCurrentType());
            std::ostringstream oss;
            oss << "Unexpected " << HTS_TYPE::label(htsType) << " order:\n"
                << "\tInput-record with pos/type/index: "
                << (getCurrentPos()+1) << "/" << HTS_TYPE::label(htsType) << "/" << getCurrentIndex()
                << " follows pos/type/index: "
                << (last.pos+1) << "/" << HTS_TYPE::label(getHtsType(last.order)) << "/" << getUserIndex(last.order) << "";
            BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
        }
    }
    else
    {
        _isStreamBegin = true;
    }

    if (_streamQueue.empty())
    {
        _isStreamEnd = true;
    }

    return (! _isStreamEnd);
}


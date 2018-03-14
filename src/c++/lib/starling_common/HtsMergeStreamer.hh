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

#pragma once

#include "blt_util/blt_types.hh"
#include "htsapi/bam_streamer.hh"
#include "htsapi/bed_streamer.hh"
#include "htsapi/vcf_streamer.hh"

#include <memory>
#include <queue>
#include <string>
#include <vector>


namespace HTS_TYPE
{
enum index_t
{
    NONE,
    BAM,
    VCF,
    BED
};

inline
const char*
label(const index_t i)
{
    switch (i)
    {
    case NONE :
        return "NONE";
    case BAM :
        return "BAM/CRAM";
    case VCF :
        return "VCF";
    case BED :
        return "BED";
    default :
        assert(false && "Unrecognized hts file type");
        return nullptr;
    }
}

template <typename T> index_t getStreamType();
template <> inline index_t getStreamType<bam_streamer>()
{
    return BAM;
}
template <> inline index_t getStreamType<vcf_streamer>()
{
    return VCF;
}
template <> inline index_t getStreamType<bed_streamer>()
{
    return BED;
}

template <typename T>
T* htsTypeFactory(
    const char* htsFilename, const char* /*referenceFilename*/, const char* region, const bool /*requireNormalized*/)
{
    return new T(htsFilename, region);
}
template <>
inline bam_streamer* htsTypeFactory(
    const char* htsFilename, const char* referenceFilename, const char* region, const bool /*requireNormalized*/)
{
    return new bam_streamer(htsFilename, referenceFilename, region);
}
template <>
inline bed_streamer* htsTypeFactory(
    const char* htsFilename, const char* /*referenceFilename*/, const char* region, const bool requireNonZeroRegionLength)
{
    return new bed_streamer(htsFilename, region, requireNonZeroRegionLength);
}
template <>
inline vcf_streamer* htsTypeFactory(
    const char* htsFilename, const char* /*referenceFilename*/, const char* region, const bool requireNormalized)
{
    return new vcf_streamer(htsFilename, region, requireNormalized);
}
}


/// An object which can register various htslib files and stream the merged output of all registered files for
/// a specified genomic region.
struct HtsMergeStreamer
{
    explicit
    HtsMergeStreamer(
        const std::string& referenceFilename);

    /// register* methods:
    ///
    /// used to register an hts data source to the streamer
    ///
    /// can be called before the first call to next(), after next() is called, these methods will
    /// generate a runtime error
    ///
    /// input index will be reflected back out by the MergeStreamer, it does not need to be unique
    ///
    /// registration order will be used to order all inputs with the same position
    ///
    const bam_streamer&
    registerBam(
        const std::string& bamFilename,
        const unsigned index = 0)
    {
        return registerHtsType(bamFilename, index,_data._bam);
    }

    /// \param[in] requireNonZeroRegionLength If true, an exception is thrown for any input bed record which with region
    ///                                       size of 0 or less, otherwise such records are skipped without error.
    const bed_streamer&
    registerBed(
        const std::string& bedFilename,
        const unsigned index = 0,
        const bool requireNonZeroRegionLength = true)
    {
        return registerHtsType(bedFilename,index,_data._bed, requireNonZeroRegionLength);
    }

    /// \param[in] requireNormalized if true an exception is thrown for any input variant records which are not
    ///                              normalized (see \ref vcf_streamer for definition)
    const vcf_streamer&
    registerVcf(
        const std::string& vcfFilename,
        const unsigned index = 0,
        const bool requireNormalized = true)
    {
        return registerHtsType(vcfFilename,index,_data._vcf, requireNormalized);
    }

    /// Reset the region over which all files scanned
    ///
    /// \param[in] region samtools-formatted genomic region string
    void
    resetRegion(const std::string& region);


    /// Advances to the next HTS record in the merged stream
    ///
    /// returns false if no more records exist
    ///
    /// This needs to be called once before any data is accessible. The first call to next will not discard any records.
    ///
    bool
    next();

    HTS_TYPE::index_t
    getCurrentType() const
    {
        return getHtsType(getCurrent().order);
    }

    unsigned
    getCurrentIndex() const
    {
        return getUserIndex(getCurrent().order);
    }

    pos_t
    getCurrentPos() const
    {
        return getCurrent().pos;
    }

    const bam_record&
    getCurrentBam() const
    {
        return *(getCurrentBamStreamer().get_record_ptr());
    }

    const bed_record&
    getCurrentBed() const
    {
        return *(getHtsStreamer(getCurrent().order, _data._bed).get_record_ptr());
    }

    const vcf_record&
    getCurrentVcf() const
    {
        return *(getCurrentVcfStreamer().get_record_ptr());
    }

    const bam_streamer&
    getCurrentBamStreamer() const
    {
        return getHtsStreamer(getCurrent().order, _data._bam);
    }

    const vcf_streamer&
    getCurrentVcfStreamer() const
    {
        return getHtsStreamer(getCurrent().order, _data._vcf);
    }

private:

    struct HtsData
    {
        std::vector<std::unique_ptr<bam_streamer>> _bam;
        std::vector<std::unique_ptr<vcf_streamer>> _vcf;
        std::vector<std::unique_ptr<bed_streamer>> _bed;
    };

    /// this sets the sort order for different HTS stream types
    /// and instances
    ///
    struct HtsRecordSortData
    {
        HtsRecordSortData(
            const pos_t initPos = 0,
            const unsigned initOrder = 0)
            :  pos(initPos), order(initOrder)
        {}

        // reverse logic implied by operator< such that the 'lower' values
        // we'd like to see first will come up on top of the
        // priority_queue
        //
        bool
        operator<(const HtsRecordSortData& rhs) const
        {
            if (pos > rhs.pos) return true;
            if (pos != rhs.pos) return false;
            return (order > rhs.order);
        }

        pos_t pos;
        // record the submission order:
        unsigned order;
    };

    struct OrderData
    {
        OrderData(
            const HTS_TYPE::index_t initHtsType,
            const unsigned initUserIndex,
            const unsigned initHtsTypeIndex)
            : htsType(initHtsType)
            , userIndex(initUserIndex)
            , htsTypeIndex(initHtsTypeIndex)
        {}

        HTS_TYPE::index_t htsType = HTS_TYPE::NONE;
        unsigned userIndex = 0;
        unsigned htsTypeIndex = 0;
    };

    const char*
    getRegionPtr()
    {
        if (_region.empty()) return nullptr;
        return _region.c_str();
    }

    const HtsRecordSortData&
    getCurrent() const
    {
        assert(_isStreamBegin && (! _isStreamEnd));
        return _streamQueue.top();
    }

    HTS_TYPE::index_t
    getHtsType(const unsigned orderIndex) const
    {
        return _order[orderIndex].htsType;
    }

    unsigned
    getUserIndex(const unsigned orderIndex) const
    {
        return _order[orderIndex].userIndex;
    }

    /// \param[in] isHighStringencyMode If true an exception is thrown for any input records which do not meet a high
    ///                                 stringency validation criteria. This criteria is different for each HTS record type.
    template <typename T>
    const T&
    registerHtsType(
        const std::string& htsFilename,
        const unsigned index,
        std::vector<std::unique_ptr<T>>& htsStreamerVec,
        const bool isHighStringencyMode = false)
    {
        static const HTS_TYPE::index_t htsType(HTS_TYPE::getStreamType<T>());
        assert(! _isStreamBegin);
        const unsigned htsTypeIndex(htsStreamerVec.size());
        const unsigned orderIndex(_order.size());
        htsStreamerVec.emplace_back(HTS_TYPE::htsTypeFactory<T>(htsFilename.c_str(), _referenceFilename.c_str(), getRegionPtr(), isHighStringencyMode));
        _order.emplace_back(htsType, index, htsTypeIndex);
        queueItem(orderIndex);
        return *(htsStreamerVec.back());
    }

    const OrderData&
    getOrderForType(
        const unsigned orderIndex,
        const HTS_TYPE::index_t expectedHtsType) const;

    template <typename T>
    T&
    getHtsStreamer(
        const unsigned orderIndex,
        std::vector<std::unique_ptr<T>>& htsStreamerVec)
    {
        static const HTS_TYPE::index_t htsType(HTS_TYPE::getStreamType<T>());
        const OrderData& orderData(getOrderForType(orderIndex, htsType));
        return *(htsStreamerVec[orderData.htsTypeIndex]);
    }

    template <typename T>
    const T&
    getHtsStreamer(
        const unsigned orderIndex,
        const std::vector<std::unique_ptr<T>>& htsStreamerVec) const
    {
        static const HTS_TYPE::index_t htsType(HTS_TYPE::getStreamType<T>());
        const OrderData& orderData(getOrderForType(orderIndex, htsType));
        return *(htsStreamerVec[orderData.htsTypeIndex]);
    }

    void
    queueItem(const unsigned orderIndex);

    /// attempt to queue an item from the same order
    /// as the head of the queue:
    void
    queueNextItem()
    {
        queueItem(getCurrent().order);
    }


    /////// data:
    std::string _referenceFilename;
    std::string _region;
    HtsData _data;
    std::vector<OrderData> _order;

    bool _isStreamBegin = false;
    bool _isStreamEnd = false;
    typedef std::priority_queue<HtsRecordSortData> queue_t;
    queue_t _streamQueue;
};

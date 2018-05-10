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

#include "appstats/RunStatsManager.hh"
#include "blt_common/map_level.hh"
#include "blt_util/depth_buffer.hh"
#include "blt_util/depth_stream_stat_range.hh"
#include "blt_util/pos_processor_base.hh"
#include "blt_util/RegionTracker.hh"
#include "blt_util/stage_manager.hh"
#include "blt_util/window_util.hh"
#include "starling_common/indel_set.hh"
#include "starling_common/IndelBuffer.hh"
#include "starling_common/PileupCleaner.hh"
#include "starling_common/pos_basecall_buffer.hh"
#include "starling_common/read_mismatch_info.hh"
#include "starling_common/starling_base_shared.hh"
#include "LocalRegionStats.hh"
#include "starling_common/starling_read_buffer.hh"
#include "starling_common/starling_streams_base.hh"
#include "starling_common/ActiveRegionDetector.hh"


#include "boost/utility.hpp"

#include <memory>
#include <string>

struct diploid_genotype;


/// \brief accumulate sequential position specific information and
/// send to a snp-calling/indel-calling routine after all position
/// information is found
///
/// pos_processor assumes that information related to each position
/// will be available in an approximately sequential fashion, where all
/// position values submitted after position X will be greater than
/// X-POS_BUFFER_SIZE+1. A violation of this assumption will trigger a
/// runtime error.
///
/// The implementation should be split into a generic semi-sequential
/// position-processor and an object with the application-specific
/// position code (here, snp-calling). Maybe once a second application
/// comes up...
///
///
/// range policy:
///
/// if begin_pos is not specified, then event processiing and
/// reporting start at the first pos >= 0 with snp/indel information
/// submitted, else at begin_pos
///
/// if end_pos is not specified, then event processing and reporting
/// end after last_pos with snp/indel information submitted, else at
/// end_pos. If the reference sequence is set then end_pos must be <
/// ref_size and will be adjusted accordingly.
///
/// Submission of snps < first_pos will be ignored. Submission of indels
/// between first_pos and first_pos-MAX_READ_SIZE will be processed but
/// not reported.
///
struct starling_pos_processor_base : public pos_processor_base, private boost::noncopyable
{
    typedef pos_processor_base base_t;

    /// \param[in] statsManager used to track aggregate pos/read/variant stats
    starling_pos_processor_base(
        const starling_base_options& opt,
        const starling_base_deriv_options& dopt,
        const reference_contig_segment& ref,
        const starling_streams_base& fileStreams,
        const unsigned sampleCount,
        RunStatsManager& statsManager);

    ~starling_pos_processor_base() override;

    /// finish position report and reset structure to ground state:
    ///
    virtual void reset();

    /// note that indel position should be normalized before calling:
    ///
    void
    insert_indel(const IndelObservation& obs,
                 const unsigned sampleIndex);

    bool is_active_region_detector_enabled()
    {
        return _opt.isHaplotypingEnabled;
    }

    /// \return reference to the read buffer from the active region detector
    ActiveRegionReadBuffer&
    getActiveRegionReadBuffer(const unsigned sampleIndex)
    {
        return _getActiveRegionDetector().getReadBuffer(sampleIndex);
    }

    /// in range [begin,end), is the estimated depth always below
    /// depth?
    bool
    is_estimated_depth_range_ge_than(
        const pos_t begin,
        const pos_t end,
        const unsigned depth,
        const unsigned sample_no) const;

    /// insert read into read buffer
    ///
    /// \return true if the alignment is accepted into the buffer (alignments can fail a number of quality checks --
    /// such as being located too far away from other alignments of the same read or having an indel that is too large.
    /// If true, the return value provides the read's id in this structure's ead buffer
    ///
    boost::optional<align_id_t>
    insert_read(
        const bam_record& br,
        const alignment& al,
        const char* chrom_name,
        const MAPLEVEL::index_t maplev,
        const unsigned sampleIndex);

    /// snv gt and stats must be reported for this pos (note only honored in strelka and starling right now)
    void
    insert_forced_output_pos(const pos_t pos);

    /// specify ploidy of region for one sample (only 0 or 1 are used for ploidy now)
    /// \returns false if this conflicts with an existing region
    bool
    insert_ploidy_region(
        const unsigned sampleIndex,
        const known_pos_range2& range,
        const unsigned ploidy);

    /// specify a region as eligible for variant calling output
    void
    insertCallRegion(
        const known_pos_range2& range);

    void
    set_head_pos(const pos_t pos);

protected:
    /// Reset current report region -- must be called before inserting region data
    ///
    /// Note that pos_processor classes only take a const ref to the reference sequence, and do not reset
    /// the reference even though this needs to be reset to the new region somewhere.
    void
    resetRegionBase(
        const std::string& chromName,
        const known_pos_range2& reportRange);

    bool
    isChromNameInitialized() const
    {
        return (! _chromName.empty());
    }

    struct LocalRegionStatsCollection
    {
        LocalRegionStatsCollection()
            : _max_winsize(0)
            , _is_last_pos(false)
            , _last_insert_pos(false)
        {}

        /// return index of new window:
        unsigned
        addNewLocalRegionStatsSize(
            const unsigned winsize)
        {
            _localRegionStatsCollection.emplace_back(winsize);
            if (winsize>_max_winsize) _max_winsize=winsize;
            return (_localRegionStatsCollection.size()-1);
        }

        /// reset average data for each window, but leave the window/winsize info in place
        void
        resetRegion()
        {
            _is_last_pos = false;
            _last_insert_pos = false;
            for (auto& localRegionStats : _localRegionStatsCollection)
            {
                localRegionStats.reset();
            }
        }

        void
        insert(const pos_t pos,
               const unsigned usedBasecallCount,
               const unsigned unusedBasecallCount,
               const unsigned spanningDeletionReadCount,
               const unsigned submappedReadCount)
        {
            if (_localRegionStatsCollection.empty()) return;
            check_skipped_pos(pos);
            insert_impl(usedBasecallCount,unusedBasecallCount,spanningDeletionReadCount,submappedReadCount);
        }

        // TODO: does check_skipped_pos need to be called here as well?
        void
        insert_null(const pos_t pos)
        {
            check_skipped_pos(pos);
            for (auto& localRegionStats : _localRegionStatsCollection)
            {
                localRegionStats.regionUsedBasecallCount.insert_null();
                localRegionStats.regionUnusedBasecallCount.insert_null();
                localRegionStats.regionSpanningDeletionReadCount.insert_null();
                localRegionStats.regionSubmappedReadCount.insert_null();
            }
        }

        const LocalRegionStats&
        getLocalRegionStats(const unsigned regionStatsIndex) const
        {
            return _localRegionStatsCollection[regionStatsIndex];
        }

    private:

        void
        check_skipped_pos(const pos_t pos)
        {
            if (_is_last_pos && (pos>(_last_insert_pos+1)))
            {
                const unsigned rep(std::min(static_cast<pos_t>(_max_winsize),(pos-(_last_insert_pos+1))));
                for (unsigned i(0); i<rep; ++i) insert_impl(0,0,0,0);
            }
            _last_insert_pos=pos;
            _is_last_pos=true;
        }

        void
        insert_impl(const unsigned usedBasecallCount,
                    const unsigned unusedBasecallCount,
                    const unsigned spanningDeletionReadCount,
                    const unsigned submappedReadCount)
        {
            for (auto& localRegionStats : _localRegionStatsCollection)
            {
                localRegionStats.regionUsedBasecallCount.insert(usedBasecallCount);
                localRegionStats.regionUnusedBasecallCount.insert(unusedBasecallCount);
                localRegionStats.regionSpanningDeletionReadCount.insert(spanningDeletionReadCount);
                localRegionStats.regionSubmappedReadCount.insert(submappedReadCount);
            }
        }

        /// Holds multiple sets of LocalRegionStats, possibly covering multiple regions sizes
        std::vector<LocalRegionStats> _localRegionStatsCollection;
        unsigned _max_winsize;
        bool _is_last_pos;
        pos_t _last_insert_pos;
    };


public:
    /// Consolidates all sample-specific data
    struct sample_info
    {
        sample_info(
            const starling_base_options& opt,
            const reference_contig_segment& ref,
            read_id_counter* ricp)
            : basecallBuffer(ref)
            , readBuffer(ricp)
            , sampleOptions(opt)
            , localRegionStatsCollection()
        {}

        void
        resetRegion()
        {
            basecallBuffer.clear();
            readBuffer.clear();
            estdepth_buff.clear();
            estdepth_buff_tier2.clear();
            localRegionStatsCollection.resetRegion();
            cleanedPileup.clear();
            ploidyRegions.clear();
        }

        pos_basecall_buffer basecallBuffer;
        starling_read_buffer readBuffer;

        /// An early estimate of read depth before realignment
        depth_buffer estdepth_buff;

        /// An early estimate of tier2 read depth before realignment
        depth_buffer estdepth_buff_tier2;

        starling_sample_options sampleOptions;

        /// Collection of local region statistics
        LocalRegionStatsCollection localRegionStatsCollection;

        // keep a single copy of this struct to reuse for every site to lower alloc costs:
        CleanedPileup cleanedPileup;

        /// track expected ploidy within this sample:
        RegionPayloadTracker<unsigned> ploidyRegions;
    };

    sample_info&
    sample(const unsigned sampleIndex = 0)
    {
        assert(sampleIndex < getSampleCount());
        return *_sample[sampleIndex];
    }

    const sample_info&
    sample(const unsigned sampleIndex = 0) const
    {
        assert(sampleIndex < getSampleCount());
        return *_sample[sampleIndex];
    }

protected:
    unsigned
    getSampleCount() const
    {
        return _sample.size();
    }

    IndelBuffer&
    getIndelBuffer()
    {
        return _indelBuffer;
    }

    const IndelBuffer&
    getIndelBuffer() const
    {
        return _indelBuffer;
    }

    CandidateSnvBuffer&
    getCandidateSnvBuffer()
    {
        return _candidateSnvBuffer;
    }

    const CandidateSnvBuffer&
    getCandidateSnvBuffer() const
    {
        return _candidateSnvBuffer;
    }

public:
    ///////////////////////////////
    // static methods:
    //
    static
    void
    report_stream_stat(const depth_stream_stat_range& ss,
                       const char* label,
                       const pos_range& pr,
                       std::ostream& os);

protected:
    bool
    is_pos_reportable(const pos_t pos) const
    {
        return _reportRange.is_pos_intersect(pos);
    }

    bool
    is_pos_preceding_reportable_range(const pos_t pos) const
    {
        return (pos < _reportRange.begin_pos());
    }

    const ActiveRegionDetector&
    getActiveRegionDetector() const
    {
        assert (_activeRegionDetector);
        return *_activeRegionDetector;
    }

private:
    ActiveRegionDetector&
    _getActiveRegionDetector() const
    {
        assert (_activeRegionDetector);
        return *_activeRegionDetector;
    }
    void
    insert_pos_submap_count(const pos_t pos,
                            const unsigned sample_no);

    void
    insert_pos_spandel_count(const pos_t pos,
                             const unsigned sample_no);

    /// update (basecall) data used to compute germline variant scoring
    /// features
    ///
    /// \param distanceFromReadEdge distance of base call from the closest read edge
    ///
    void
    updateGermlineScoringMetrics(
        const pos_t pos,
        const unsigned sample_no,
        const uint8_t call_id,
        const uint8_t qscore,
        const uint8_t mapq,
        const unsigned cycle,
        const unsigned distanceFromReadEdge,
        const bool is_submapped);

    void
    updateSomaticScoringMetrics(
        const pos_t pos,
        const unsigned sample_no,
        const bool is_teir1,
        const uint8_t call_id,
        const bool is_call_filter,
        const uint16_t readPos,
        const uint16_t readLength);

    void
    insert_pos_basecall(const pos_t pos,
                        const unsigned sample_no,
                        const bool is_tier1,
                        const base_call& bc);

    void
    process_pos(const int stage_no,
                const pos_t pos) override;

    void
    load_read_in_depth_buffer(const read_segment& rseg,
                              const unsigned sample_no);

    /// process read segment alignment for indels
    void
    init_read_segment(
        const read_segment& rseg,
        const unsigned sampleIndex);

    /// Initialize all spliced read segments buffered at the given position
    ///
    void
    initializeSplicedReadSegmentsAtPos(const pos_t pos);

    /// For all reads buffered at the current position:
    /// 1) determine the set of candidate indels that the read overlaps
    /// 2) determine the set of private indels within the read's discovery alignments
    /// 3) Find the most likely alignment given both sets of indels
    /// 4) evaluate the probability that the read supports each candidate indel vs. the reference
    ///
    void
    align_pos(const pos_t pos);

    /// adjust read buffer position so that reads are buffered in sorted
    /// order after realignment:
    ///
    void
    rebuffer_pos_reads(const pos_t pos);

    /// Add reads buffered at position into a basecall pileup
    /// to allow for downstream depth and site genotyping calculations
    ///
    void
    pileup_pos_reads(const pos_t pos);

    /// Add a single read segment into the the basecall pileup
    void
    pileup_read_segment(
        const read_segment& rseg,
        const unsigned sampleIndex);

    void
    write_reads(const pos_t pos);

    /// maintain stats for depth, etc...
    void
    process_pos_stats(
        const pos_t pos);

    /// maintain per-sample stats for depth, etc...
    void
    process_pos_sample_stats(
        const pos_t pos,
        const unsigned sample_no);

    void
    process_pos_variants(
        const pos_t pos,
        const bool isPosPrecedingReportableRange);

    //////
    virtual
    void
    process_pos_variants_impl(
        const pos_t pos,
        const bool isPosPrecedingReportableRange) = 0;

    virtual
    void
    clear_pos_annotation(const pos_t /*pos*/) {}

    virtual
    bool
    is_save_pileup_buffer()
    {
        return false;
    }

    /// (re-)initialize active region detector (this can safely be done anytime after _sample vector is initialized)
    void
    resetActiveRegionDetector();

protected:
    virtual
    void
    run_post_call_step(
        const int /*stage_no*/,
        const pos_t /*pos*/)
    {}

    unsigned
    get_largest_read_size() const
    {
        return _rmi.size();
    }

private:
    // return false if read is too large
    bool
    update_largest_read_size(const unsigned rs);

protected:
    unsigned
    get_largest_total_indel_ref_span_per_read() const
    {
        return _largest_total_indel_ref_span_per_read;
    }

private:
    void
    update_largest_indel_ref_span(const unsigned is);

    void
    update_largest_total_indel_ref_span_per_read(const unsigned is);

    void
    update_stageman();

    virtual
    void
    post_align_clear_pos(const pos_t /*pos*/) {}

    bool
    empty() const
    {
        if (! _is_skip_process_pos)
        {
            const unsigned sampleCount(getSampleCount());
            for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
            {
                const sample_info& sif(sample(sampleIndex));
                if (! sif.readBuffer.empty()) return false;
                if (! sif.basecallBuffer.empty()) return false;
            }
            if (! _indelBuffer.empty()) return false;
            if (! _candidateSnvBuffer.empty()) return false;
            if (! _forced_output_pos.empty()) return false;
            if (! derived_empty()) return false;
            _is_skip_process_pos=true;
        }
        return true;
    }

    /// allow a derived class to declare non-empty status:
    virtual
    bool
    derived_empty() const
    {
        return true;
    }

protected:

    bool
    is_forced_output_pos(const pos_t pos) const
    {
        return (_forced_output_pos.find(pos) != _forced_output_pos.end());
    }

    unsigned
    get_ploidy(
        const pos_t pos,
        const unsigned sampleIndex) const
    {
        const auto val(sample(sampleIndex).ploidyRegions.isIntersectRegion(pos));
        return (val ? *val : 2u);
    }

    //////////////////////////////////
    // data:
    //
    const starling_base_options& _opt;
    const starling_base_deriv_options& _dopt;
    const reference_contig_segment& _ref;
    const starling_streams_base& _streams;
    RunStatsManager& _statsManager;

    // read-length data structure used to compute mismatch density filter:
    read_mismatch_info _rmi;

    /// Largest delete length observed for any one indel (but not greater than max_delete_size)
    unsigned _largest_indel_ref_span;

    // largest
    unsigned _largest_total_indel_ref_span_per_read;

    std::unique_ptr<stage_manager> _stagemanPtr;

    std::string _chromName;
    known_pos_range2 _reportRange;

    /// when callRegions mode is on, calling is restricted to the regions defined in this object
    RegionTracker _callRegions;

    /// used to keep read id's unique across multiple samples:
    read_id_counter _ric;

    std::vector<std::unique_ptr<sample_info>> _sample;

    std::unique_ptr<diploid_genotype> _empty_dgt[N_BASE];

    std::set<pos_t> _forced_output_pos;

    PileupCleaner _pileupCleaner;

private:
    IndelBuffer _indelBuffer;
    CandidateSnvBuffer _candidateSnvBuffer;
    std::unique_ptr<ActiveRegionDetector> _activeRegionDetector;
};

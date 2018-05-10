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

#ifdef _MSC_VER
#pragma warning(disable:4355)
#endif

#include "starling_pos_processor_indel_util.hh"
#include "starling_read_align.hh"
#include "starling_read_util.hh"

#include "blt_common/position_snp_call_pprob_digt.hh"
#include "blt_util/depth_buffer_util.hh"
#include "blt_util/io_util.hh"
#include "blt_util/log.hh"
#include "htsapi/bam_seq_read_util.hh"
#include "starling_common/AlleleReportInfo.hh"

#include <iomanip>


//#define DEBUG_PPOS

// largest_read_size grows dynamically with observed read size, this
// is used to initialize the setting prior to observing any reads:
//
// initial setting is large to help consistently deal with grouperisms:
//
const unsigned STARLING_INIT_LARGEST_READ_SIZE(250);
const double STARLING_LARGEST_READ_SIZE_PAD(1.25);

// additional distance between HEAD and READ_BUFFER for haplotyping
const unsigned HAPLOTYPING_PADDING(200);

// largest indel_size grows dynamically with observed indel size until
// hitting maxIndelSize. Initialized to the follow value prior to
// observation:



//////////////////////////////////////////////
// file-specific static functions:
//



static
void
report_pos_range(const pos_range& pr,
                 std::ostream& os)
{
    // convert pos_range to 1-indexed inclusive interval for output:
    os << "begin: ";
    if (pr.is_begin_pos)
    {
        os << pr.begin_pos+1;
    }
    else
    {
        os << "NA";
    }

    os << " end: ";
    if (pr.is_end_pos)
    {
        os << pr.end_pos;
    }
    else
    {
        os << "NA";
    }
}



static
unsigned
get_read_buffer_size(const unsigned largest_read_size,
                     const unsigned largest_total_indel_span_per_read)
{
    return (largest_read_size+largest_total_indel_span_per_read);
}



// static methods:
//
void
starling_pos_processor_base::
report_stream_stat(const depth_stream_stat_range& ss,
                   const char* label,
                   const pos_range& pr,
                   std::ostream& os)
{
    os << label << " ";
    report_pos_range(pr,os);
    os << " " << ss << "\n";
}




////////////////////////////////////////////////////////
// define starling_pos_processor_base stages:
//

namespace STAGE
{

// stage into which pileup entries must fit:
static
int
get_pileup_stage_no(const starling_base_options& /*opt*/)
{
    return (static_cast<int>(POST_ALIGN));
}

// given largest read, and indel ref span per read, get stage
// lengths:
//
static
stage_data
get_stage_data(
    const unsigned largest_read_size,
    const unsigned largest_total_indel_ref_span_per_read,
    const starling_base_options& opt,
    const starling_base_deriv_options& dopt)
{
    stage_data sdata;

    // HEAD contains everything before the head position in the
    // stage processing pipeline, this should be:
    //
    // 1) Exon entries that extend past the head position
    //
    // HEAD has no defined size -- it is everything stored before
    // the head position
    //
    sdata.add_stage(HEAD);


    // READ_BUFFER is where most normal genomic reads are read in
    // and processed for indels, it needs enough room to collect the full
    // candidate indel map before we start read alignment:
    //
    // occurring at the beginning of stage READ_BUFFER (or before):
    // collect and buffer read
    // enter read into estimated depth
    //
    // occurring at the end of stage READ_BUFFER:
    // read realignment
    // base-call pos entry (pileup)
    //
    sdata.add_stage(READ_BUFFER,HEAD,get_read_buffer_size(largest_read_size,largest_total_indel_ref_span_per_read)+HAPLOTYPING_PADDING);

    //
    // POST_ALIGN_STAGE - the end of this stage is where snp and indel
    // calling occur. It needs to be separated from the end of
    // READ_BUFFER_STAGE with enough room to allow for the repositioning
    // of reads to the 'left' during re-alignment.
    //
    // At the end of this stage we conduct:
    //
    // indel calling
    // snp-calling
    // realigned read output
    //
    sdata.add_stage(POST_ALIGN,READ_BUFFER,largest_total_indel_ref_span_per_read);

    sdata.add_stage(CLEAR_SITE_ANNOTATION,POST_ALIGN,largest_total_indel_ref_span_per_read);

    sdata.add_stage(CLEAR_READ_BUFFER,POST_ALIGN,0);

    // CLEAR_INDEL_BUFFER is the point at which all indel scoring data is cleared
    //
    // it needs to persistent long enough after indel scoring so that overlapping indels in
    // 5' of the indel/allele currently being called are still available.
    //
    sdata.add_stage(CLEAR_INDEL_BUFFER,POST_ALIGN,largest_total_indel_ref_span_per_read);

    // POST_CALL is used to free up the pileup data. This data is
    // preserved for one extra cycle after snp and indel calling so that
    // the indel caller can look one base 'to the left' of its call
    // location (excepting right breakpoints) to get the depth for the
    // current indel call
    //
    // At the end of this stage we free the position's pileup buffer
    //
    static const unsigned pileup_free_delay(1);
    sdata.add_stage(POST_CALL,POST_ALIGN,pileup_free_delay);

    // dynamic stages after POST_CALL are used to region based information around a site
    //
    const int pileupStageIndex(get_pileup_stage_no(opt));

    const auto& postCallStages(dopt.getPostCallStage());
    const unsigned postCallStageCount(postCallStages.size());
    for (unsigned postCallStageIndex(0); postCallStageIndex<postCallStageCount; ++postCallStageIndex)
    {
        sdata.add_stage(SIZE+postCallStageIndex,pileupStageIndex,postCallStages[postCallStageIndex]);
    }

    return sdata;
}
}

starling_pos_processor_base::
starling_pos_processor_base(
    const starling_base_options& opt,
    const starling_base_deriv_options& dopt,
    const reference_contig_segment& ref,
    const starling_streams_base& fileStreams,
    const unsigned sampleCount,
    RunStatsManager& statsManager)
    : base_t()
    , _opt(opt)
    , _dopt(dopt)
    , _ref(ref)
    , _streams(fileStreams)
    , _statsManager(statsManager)
    , _rmi(STARLING_INIT_LARGEST_READ_SIZE)
    , _largest_indel_ref_span(opt.maxIndelSize)
    , _largest_total_indel_ref_span_per_read(_largest_indel_ref_span)
    , _sample(sampleCount)
    , _pileupCleaner(opt)
    , _indelBuffer(opt,dopt,ref)
    , _candidateSnvBuffer(sampleCount)
{
    assert(sampleCount != 0);

    for (auto& sampleVal : _sample)
    {
        sampleVal.reset(new sample_info(_opt, ref, &_ric));
    }

    // this can safely be called after initializing _sample above
    resetActiveRegionDetector();

    if (_opt.is_all_sites())
    {
        // pre-calculate qscores for sites with no observations:
        //
        snp_pos_info good_pi;
        static const std::vector<float> dependent_eprob;
        const extended_pos_info good_epi(good_pi,dependent_eprob);
        for (unsigned b(0); b<N_BASE; ++b)
        {
            good_pi.set_ref_base(id_to_base(b));
            _empty_dgt[b].reset(new diploid_genotype);
            _dopt.pdcaller().position_snp_call_pprob_digt(_opt,good_epi,
                                                          *_empty_dgt[b],
                                                          _opt.is_all_sites());
        }
    }
}


void
starling_pos_processor_base::
resetActiveRegionDetector()
{
    _activeRegionDetector.reset(
        new ActiveRegionDetector(
            _ref, _indelBuffer, _candidateSnvBuffer,
            _opt.maxIndelSize, getSampleCount(), _opt.isSomaticCallingMode));
}



void
starling_pos_processor_base::
update_stageman()
{
    _stagemanPtr->revise_stage_data(
        STAGE::get_stage_data(get_largest_read_size(),
                              get_largest_total_indel_ref_span_per_read(),
                              _opt,
                              _dopt));
}


bool
starling_pos_processor_base::
update_largest_read_size(const unsigned rs)
{
    if (rs>STRELKA_MAX_READ_SIZE) return false;

    if (rs<=get_largest_read_size()) return true;
    _rmi.resize(rs);
    update_stageman();
    return true;
}



void
starling_pos_processor_base::
update_largest_indel_ref_span(const unsigned is)
{
    if (is<=_largest_indel_ref_span) return;
    assert(is<=_opt.maxIndelSize);
    _largest_indel_ref_span=std::min(is,_opt.maxIndelSize);
    update_largest_total_indel_ref_span_per_read(is);
    update_stageman();
}



void
starling_pos_processor_base::
update_largest_total_indel_ref_span_per_read(const unsigned is)
{
    if (is<=_largest_total_indel_ref_span_per_read) return;
    _largest_total_indel_ref_span_per_read=is;
    update_stageman();
}



starling_pos_processor_base::
~starling_pos_processor_base()
{
}



void
starling_pos_processor_base::
reset()
{
    if (_stagemanPtr)
    {
        _stagemanPtr->reset();
    }
    _activeRegionDetector->clear();
}



void
starling_pos_processor_base::
resetRegionBase(
    const std::string& chromName,
    const known_pos_range2& reportRange)
{
    reset();
    _chromName = chromName;
    _reportRange = reportRange;

    // note that reseting these 'largest indel seen' values shouldn't really be necessary/important,
    // but it contributes to easier verification that a series of regions put into the caller will
    // give the same result as those regions processed one at a time
    _largest_indel_ref_span = _opt.maxIndelSize;
    _largest_total_indel_ref_span_per_read = _largest_indel_ref_span;

    const pos_range pr(_reportRange.begin_pos(), _reportRange.end_pos());
    _stagemanPtr.reset(new stage_manager(STAGE::get_stage_data(STARLING_INIT_LARGEST_READ_SIZE, get_largest_total_indel_ref_span_per_read(), _opt, _dopt), pr, *this));

    _forced_output_pos.clear();
    _indelBuffer.clearIndels();
    _candidateSnvBuffer.clearSnvs();

    /// TODO, it might be better to have some kind of regionReset() on this structure
    ///  -- not clear how to do this accurately, so for now we just nuke and replace the entire object
    resetActiveRegionDetector();

    for (auto& sampleVal : _sample)
    {
        sampleVal->resetRegion();
    }
}



void
starling_pos_processor_base::
insert_indel(
    const IndelObservation& obs,
    const unsigned sampleIndex)
{
    // ppr advance is controlled by the start positions of reads and
    // relatively cheap to store (so long as we aren't including
    // gigantic insert sequences) and do not scale up linearly with
    // increased coverage like reads do. For this reason our strategy
    // is to buffer the indels as far ahead as possible while leaving
    // the read buffer sizes fixed at a smaller value.
    //

    // address STARKA-248
    //
    // This is designed to filter-out null operations from the cigar string -- that is, matched insert/delete segments that just reinsert the
    // reference. It is, in principal, better to filter here so that these spurious indels don't gum up the realignment and analysis downstream,
    // however the risk of instability is  high at this point. Assuming these occur rarely in the input bam, the impact of leaving these in the
    // process should be relatively low, and keep starka stable in the short term.
    //
    // Consider reactivating this filter in the longer term release path
#if 0
    if ((!obs.key.is_breakpoint()) && (obs.key.delete_length() == obs.key.insert_length()))
    {
        assert(obs.data.insert_seq.size() == obs.key.insert_length());

        std::string dseq;
        _ref.get_substring(obs.key.pos,obs.key.delete_length(),dseq);
        if (obs.data.insert_seq == dseq)
        {
            // choose to either filter or throw:
#if 1
            return true;
#else
            std::ostringstream oss;
            oss << "Read alignment contains indel allele matching reference obs: " << obs << "\n";
            throw blt_exception(oss.str().c_str());
#endif
        }
    }
#endif

    try
    {
        _stagemanPtr->validate_new_pos_value(obs.key.pos,STAGE::READ_BUFFER);

        // dynamically scale maximum indel size:
        const unsigned len(std::min(static_cast<unsigned>((obs.key.delete_length())),_opt.maxIndelSize));
        update_largest_indel_ref_span(len);


        if (is_active_region_detector_enabled())
        {
            getActiveRegionReadBuffer(sampleIndex).insertIndel(obs);
        }
        else
        {
            getIndelBuffer().addIndelObservation(sampleIndex, obs);
        }
        if (obs.data.is_forced_output) _is_skip_process_pos=false;
    }
    catch (...)
    {
        log_os << "Exception caught while attempting to insert indel: " << obs << "\n";
        throw;
    }
}



void
starling_pos_processor_base::
insert_forced_output_pos(const pos_t pos)
{
    _stagemanPtr->validate_new_pos_value(pos,STAGE::READ_BUFFER);
    _forced_output_pos.insert(pos);
    _is_skip_process_pos=false;
}



bool
starling_pos_processor_base::
insert_ploidy_region(
    const unsigned sampleIndex,
    const known_pos_range2& range,
    const unsigned ploidy)
{
    assert(ploidy==0 || ploidy==1);
    _stagemanPtr->validate_new_pos_value(range.begin_pos(),STAGE::READ_BUFFER);
    return sample(sampleIndex).ploidyRegions.addRegion(range,ploidy);
}



void
starling_pos_processor_base::
insertCallRegion(
    const known_pos_range2& range)
{
    _stagemanPtr->validate_new_pos_value(range.begin_pos(),STAGE::READ_BUFFER);
    _callRegions.addRegion(range);
    _is_skip_process_pos=false;
}



bool
starling_pos_processor_base::
is_estimated_depth_range_ge_than(
    const pos_t begin,
    const pos_t end,
    const unsigned depth,
    const unsigned sample_no) const
{
    const auto& est1(sample(sample_no).estdepth_buff);
    const auto& est2(sample(sample_no).estdepth_buff_tier2);

    assert(begin <= end);
    for (pos_t i(begin); i<=end; ++i)
    {
        if ((est1.val(i)+est2.val(i)) >= depth) return true;
    }
    return false;
}



boost::optional<align_id_t>
starling_pos_processor_base::
insert_read(
    const bam_record& br,
    const alignment& al,
    const char* chrom_name,
    const MAPLEVEL::index_t maplev,
    const unsigned sampleIndex)
{
    boost::optional<align_id_t> retval;

    if (0 != strcmp(_chromName.c_str(),chrom_name))
    {
        log_os << "ERROR: starling_pos_processor_base.insert_read(): read has unexpected sequence name: '" << chrom_name << "' expecting: '" << _chromName << "'\n"
               << "\tread_key: " << read_key(br) << "\n";
        exit(EXIT_FAILURE);
    }

    // Check at al.pos rather than buffer pos because we need to
    // insure that all candidate indels get entered by the end of HEAD
    // stage. If the read aligns back past the end of head stage it is
    // ok so long as we know it will not generate any candidate indels
    // in this region:
    //
    /// TODO this should never happen, investigate/fix turn this into an assertion
    if (! _stagemanPtr->is_new_pos_value_valid(al.pos,STAGE::HEAD))
    {
        log_os << "WARNING: skipping alignment for read outside of buffer range: " << br << "\n";
        return retval;
    }

    starling_read_buffer& rbuff(sample(sampleIndex).readBuffer);

    // check whether the read buffer has reached max capacity
    if (_opt.isMaxBufferedReads())
    {
        if (rbuff.size() >= _opt.maxBufferedReads)
        {
            return retval;
        }
    }

    // assume that pos_processor, as a container, is no longer empty...
    _is_skip_process_pos=false;

    // update read_size:
    {
        const unsigned rs(static_cast<unsigned>(STARLING_LARGEST_READ_SIZE_PAD*br.read_size()));
        if (! update_largest_read_size(rs))
        {
            std::ostringstream oss;
            oss << "Input read size: " << br.read_size() << " exceeds maximum.";
            throw blt_exception(oss.str().c_str());
        }
    }

    // insert the read:
    retval.reset(rbuff.add_read_alignment(br,al,maplev));

    // must initialize initial read_segments "by-hand":
    //
    // TODO get this streamlined into the pos-processor
    //
    if (retval)
    {
        const starling_read* sread_ptr(rbuff.get_read(*retval));
        assert(nullptr!=sread_ptr);

        // update depth-buffer for the whole read:
        load_read_in_depth_buffer(sread_ptr->get_full_segment(),sampleIndex);

        // update other data for only the first read segment
        const seg_id_t firstReadSegmentIndex(sread_ptr->isSpliced() ? 1 : 0 );
        init_read_segment(sread_ptr->get_segment(firstReadSegmentIndex),sampleIndex);
    }

    return retval;
}



static
INDEL_ALIGN_TYPE::index_t
translate_maplev_to_indel_type(const MAPLEVEL::index_t i)
{
    switch (i)
    {
    case MAPLEVEL::TIER1_MAPPED:
        return INDEL_ALIGN_TYPE::GENOME_TIER1_READ;
    case MAPLEVEL::TIER2_MAPPED:
        return INDEL_ALIGN_TYPE::GENOME_TIER2_READ;
    case MAPLEVEL::SUB_MAPPED:
        return INDEL_ALIGN_TYPE::GENOME_SUBMAP_READ;
    default:
        log_os << "ERROR: unexpected maplevel: " << MAPLEVEL::get_label(i) << "\n";
        exit(EXIT_FAILURE);
    }
}



// only acts on genomic mapped reads:
void
starling_pos_processor_base::
load_read_in_depth_buffer(const read_segment& rseg,
                          const unsigned sample_no)
{
    const alignment& al(rseg.getInputAlignment());
    if (al.empty()) return;

    const MAPLEVEL::index_t maplev(rseg.getInputAlignmentMapLevel());
    const bool is_usable_mapping(MAPLEVEL::TIER1_MAPPED == maplev);
    if (is_usable_mapping)
    {
        add_alignment_to_depth_buffer(al.pos,al.path,sample(sample_no).estdepth_buff);
    }
    else if (maplev == MAPLEVEL::TIER2_MAPPED)
    {
        add_alignment_to_depth_buffer(al.pos,al.path,sample(sample_no).estdepth_buff_tier2);
    }
}



void
starling_pos_processor_base::
init_read_segment(
    const read_segment& rseg,
    const unsigned sampleIndex)
{
    const alignment& al(rseg.getInputAlignment());
    if (al.empty()) return;

    const MAPLEVEL::index_t maplev(rseg.getInputAlignmentMapLevel());
    const INDEL_ALIGN_TYPE::index_t iat(translate_maplev_to_indel_type(maplev));

    const bam_seq bseq(rseg.get_bam_read());
    try
    {
        const unsigned total_indel_ref_span_per_read =
            addAlignmentIndelsToPosProcessor(_opt.maxIndelSize, _ref,
                                             al, bseq, *this, iat, rseg.getReadIndex(), sampleIndex, rseg.get_segment_edge_pin(),
                                             rseg.map_qual() == 0);

        update_largest_total_indel_ref_span_per_read(total_indel_ref_span_per_read);
    }
    catch (...)
    {
        log_os << "\nException caught in addAlignmentIndelsToPosProcessor() while processing record: " << rseg << "\n";
        throw;
    }
}



void
starling_pos_processor_base::
initializeSplicedReadSegmentsAtPos(const pos_t pos)
{
    const unsigned sampleCount(getSampleCount());
    for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
    {
        read_segment_iter ri(sample(sampleIndex).readBuffer.get_pos_read_segment_iter(pos));
        for (read_segment_iter::ret_val r; true; ri.next())
        {
            r=ri.get_ptr();
            if (nullptr==r.first) break;
            // full_segments of unspliced reads and the initial
            // segment of spliced reads are initialized outside of the
            // process_pos framework, so this routine only initializes
            // the second segment or higher:
            //
            if (r.second<2) continue;
            init_read_segment(r.first->get_segment(r.second),sampleIndex);
        }
    }
}



/// get max-min bounds in which reads can be realigned:
static
known_pos_range
get_realignment_range(const pos_t pos,
                      const stage_data& sdata)
{
    const unsigned head_offset(sdata.get_stage_id_shift(STAGE::HEAD));
    const unsigned buffer_offset(sdata.get_stage_id_shift(STAGE::READ_BUFFER));
    const unsigned post_offset(sdata.get_stage_id_shift(STAGE::POST_ALIGN));
    assert(buffer_offset>head_offset);
    assert(post_offset>buffer_offset);

    pos_t min_offset(static_cast<pos_t>(post_offset-buffer_offset));
    pos_t max_offset(static_cast<pos_t>(buffer_offset-head_offset));

    // shrink values by one to be safe. We aren't allowed to step outside of the realign boundary at all
    /// \TODO get exact answer on buffer range boundary so that uncertainty is eliminated here:
    min_offset = std::max(0,min_offset-1);
    max_offset = std::max(0,max_offset-1);

    const pos_t min_pos(std::max(static_cast<pos_t>(0),pos-min_offset));
    const pos_t max_pos(pos+1+max_offset); // +1 to follow right-open range convention

    return known_pos_range(min_pos, max_pos);
}



void
starling_pos_processor_base::
align_pos(const pos_t pos)
{
    known_pos_range realign_buffer_range(get_realignment_range(pos, _stagemanPtr->get_stage_data()));

    const unsigned sampleCount(getSampleCount());
    for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
    {
        sample_info& sif(sample(sampleIndex));
        read_segment_iter ri(sif.readBuffer.get_pos_read_segment_iter(pos));
        for (read_segment_iter::ret_val r; true; ri.next())
        {
            r=ri.get_ptr();
            if (nullptr == r.first) break;
            read_segment& rseg(r.first->get_segment(r.second));
            if (not (_opt.is_realign_submapped_reads || rseg.is_tier1or2_mapping())) continue;

            try
            {
                realignAndScoreRead(_opt, _dopt, sif.sampleOptions, _ref, realign_buffer_range, sampleIndex,
                                    rseg, getIndelBuffer());
            }
            catch (...)
            {
                log_os << "Exception caught in align_pos() while realigning segment: "
                       << static_cast<int>(r.second) << " of read: " << (*r.first) << "\n";
                throw;
            }
            // check that read has not been realigned too far to the left:
            if (rseg.is_realigned)
            {
                if (! _stagemanPtr->is_new_pos_value_valid(rseg.realignment.pos,STAGE::POST_ALIGN))
                {
                    log_os << "WARNING: read realigned outside bounds of realignment stage buffer. Skipping...\n"
                           << "\tread: " << rseg.key() << "\n";
                    rseg.is_invalid_realignment=true;
                }
            }
        }
    }
}



void
starling_pos_processor_base::
set_head_pos(const pos_t pos)
{
    _stagemanPtr->validate_new_pos_value(pos,STAGE::READ_BUFFER);
    _stagemanPtr->handle_new_pos_value(pos);
}



void
starling_pos_processor_base::
process_pos(const int stage_no,
            const pos_t pos)
{
#if 0
    log_os << "pos,stage_no: " << pos << " " << stage_no << "\n";
#endif

    assert(isChromNameInitialized());

    if (empty()) return;

    const unsigned sampleCount(getSampleCount());

    if        (stage_no==STAGE::HEAD)
    {
        initializeSplicedReadSegmentsAtPos(pos);
        if (is_active_region_detector_enabled())
        {
            _getActiveRegionDetector().updateEndPosition(pos);
        }
    }
    else if (stage_no==STAGE::READ_BUFFER)
    {
        align_pos(pos);
        pileup_pos_reads(pos);
        write_reads(pos);

        if (is_active_region_detector_enabled())
        {
            _getActiveRegionDetector().clearReadBuffer(pos);
        }
    }
    else if (stage_no==STAGE::POST_ALIGN)
    {
        static const pos_t precedingCallRegionDistance(200);

        // start variant calling if pos is either:
        // (1) in a call region
        // (2) no more than precedingCallRegionDistance bases upstream of a call region
        //
        // All variant callers run in case (1). In case (2) we only run variant callers which accumulate history
        // from position to position. At this time, the only such variant caller is the germline indel caller,
        // which stores history over muliptle calls in order to (more) correctly handle overlapping indels.
        //
        // Note here that "reportable range" refers the region of responsible for this variant calling process
        // but "call region" may represent several subrgions of the reportable range if a call regions BED file
        // is given (otherwise "call region" is equal to "reportable range")

        const bool isPosInPrimaryReportingRegion(is_pos_reportable(pos));

        /// is this position reportable in the VCF output?
        bool isPosReportable(false);

        /// is this position upstream of one that is reportable in the VCF output by no more than
        /// precedingCallRegionDistance bases?
        ///
        bool isPosPrecedingReportable(false);

        if (not isPosInPrimaryReportingRegion)
        {
            isPosPrecedingReportable = (is_pos_preceding_reportable_range(pos) and is_pos_reportable(pos+precedingCallRegionDistance));
            if (_opt.isUseCallRegions())
            {
                if (isPosPrecedingReportable)
                {
                    isPosPrecedingReportable = _callRegions.isIntersectRegion(known_pos_range2((pos-precedingCallRegionDistance),pos));
                }
            }
        }
        else
        {
            if (_opt.isUseCallRegions())
            {
                isPosReportable = _callRegions.isIntersectRegion(pos);
                if (not isPosReportable)
                {
                    isPosPrecedingReportable = _callRegions.isIntersectRegion(known_pos_range2((pos-precedingCallRegionDistance),pos));
                }
            }
            else
            {
                isPosReportable = true;
            }
        }

        if (isPosReportable or isPosPrecedingReportable)
        {
            process_pos_variants(pos, isPosPrecedingReportable);
        }

        if (is_active_region_detector_enabled())
        {
            _getActiveRegionDetector().clearPosToActiveRegionIdMapUpToPos(pos);
            for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
            {
                _candidateSnvBuffer.clearUpToPos(sampleIndex, pos);
            }
        }

        // everything else:
        post_align_clear_pos(pos);
    }
    else if (stage_no==STAGE::CLEAR_SITE_ANNOTATION)
    {
        _forced_output_pos.erase(pos);
        for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
        {
            sample(sampleIndex).ploidyRegions.removeToPos(pos);
        }
        clear_pos_annotation(pos);
    }
    else if (stage_no==STAGE::CLEAR_READ_BUFFER)
    {
        for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
        {
            sample(sampleIndex).readBuffer.clear_to_pos(pos);
        }
    }
    else if (stage_no==STAGE::CLEAR_INDEL_BUFFER)
    {
        getIndelBuffer().clearIndelsAtPosition(pos);
    }
    else if (stage_no==STAGE::POST_CALL)
    {
        for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
        {
            sample_info& sif(sample(sampleIndex));
            sif.estdepth_buff.clear_pos(pos);
            sif.estdepth_buff_tier2.clear_pos(pos);

            // if we are doing short-range phasing, suspend pileup buffer clear
            // while phasing block is being built:
            if (! is_save_pileup_buffer())
            {
                sif.basecallBuffer.clear_to_pos(pos);
            }
        }
    }
    else if (stage_no>=STAGE::SIZE)
    {
        run_post_call_step(stage_no, pos);
    }
    else
    {
        log_os << "ERROR: unexpected processing stage in starling_pos_processor_base\n";
        exit(EXIT_FAILURE);
    }
}



void
starling_pos_processor_base::
insert_pos_submap_count(const pos_t pos,
                        const unsigned sample_no)
{
    // assume pos has been pre-checked:

    sample(sample_no).basecallBuffer.insert_pos_submap_count(pos);
}



void
starling_pos_processor_base::
insert_pos_spandel_count(const pos_t pos,
                         const unsigned sample_no)
{
    // assume pos has been pre-checked:

    sample(sample_no).basecallBuffer.insert_pos_spandel_count(pos);
}



void
starling_pos_processor_base::
updateGermlineScoringMetrics(
    const pos_t pos,
    const unsigned sample_no,
    const uint8_t call_id,
    const uint8_t qscore,
    const uint8_t mapq,
    const unsigned cycle,
    const unsigned distanceFromReadEdge,
    const bool is_submapped)
{
    // assume pos is already valid:
    auto& bcbuff(sample(sample_no).basecallBuffer);
    bcbuff.updateGermlineScoringMetrics(_ref.get_base(pos), pos, call_id, qscore, mapq, cycle, distanceFromReadEdge, is_submapped);
}


void
starling_pos_processor_base::
updateSomaticScoringMetrics(
    const pos_t pos,
    const unsigned sample_no,
    const bool is_tier1,
    const uint8_t call_id,
    const bool is_call_filter,
    const uint16_t readPos,
    const uint16_t readLength)
{
    // assume pos is already valid:
    auto& bcbuff(sample(sample_no).basecallBuffer);
    if (is_tier1 && (sample_no != 0) && (! is_call_filter))
    {
        bcbuff.update_read_pos_ranksum(_ref.get_base(pos),pos,call_id,readPos);
        bcbuff.insert_alt_read_pos(pos,call_id,readPos,readLength);
    }
}


void
starling_pos_processor_base::
insert_pos_basecall(const pos_t pos,
                    const unsigned sample_no,
                    const bool is_tier1,
                    const base_call& bc)
{
    // assume pos is already valid:

    sample(sample_no).basecallBuffer.insert_pos_basecall(pos,is_tier1,bc);
}



#if 1
static
pos_t
get_new_read_pos(const read_segment& rseg)
{
    // get the best alignment for the read:
    const alignment* best_al_ptr(&(rseg.getInputAlignment()));
    if (rseg.is_realigned) best_al_ptr=&(rseg.realignment);

    if (best_al_ptr->empty()) return rseg.buffer_pos;     // a grouper contig read which was not realigned...
    else                     return best_al_ptr->pos;
}



void
starling_pos_processor_base::
rebuffer_pos_reads(const pos_t pos)
{
    // need to queue up read changes and run at the end so that we
    // don't invalidate read buffer iterators
    //
    typedef std::pair<std::pair<align_id_t,seg_id_t>,pos_t> read_pos_t;

    const unsigned sampleCount(getSampleCount());
    for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
    {
        sample_info& sif(sample(sampleIndex));
        std::vector<read_pos_t> new_read_pos;
        read_segment_iter ri(sif.readBuffer.get_pos_read_segment_iter(pos));
        for (read_segment_iter::ret_val r; true; ri.next())
        {
            r=ri.get_ptr();
            if (nullptr==r.first) break;
            read_segment& rseg(r.first->get_segment(r.second));

            const pos_t new_pos(get_new_read_pos(rseg));
            if ((new_pos!=pos) &&
                (_stagemanPtr->is_new_pos_value_valid(new_pos,STAGE::POST_ALIGN)))
            {
                new_read_pos.push_back(std::make_pair(std::make_pair(rseg.getReadIndex(),r.second),new_pos));
            }
        }

        const unsigned nr(new_read_pos.size());
        for (unsigned i(0); i<nr; ++i)
        {
            sif.readBuffer.rebuffer_read_segment(new_read_pos[i].first.first,
                                                 new_read_pos[i].first.second,
                                                 new_read_pos[i].second);
        }
    }
}
#endif



// this method could be const if we change the read buffer to have
// const iterators
void
starling_pos_processor_base::
write_reads(const pos_t pos)
{
    const unsigned sampleCount(getSampleCount());
    for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
    {
        bam_dumper* bamd_ptr(_streams.realign_bam_ptr(sampleIndex));
        if (nullptr == bamd_ptr) continue;
        bam_dumper& bamd(*bamd_ptr);

        read_segment_iter ri(sample(sampleIndex).readBuffer.get_pos_read_segment_iter(pos));
        read_segment_iter::ret_val r;

        while (true)
        {
            r=ri.get_ptr();
            if (nullptr==r.first) break;
            if (r.first->getExonCount()==r.second)
            {
                r.first->write_bam(bamd);
            }
            ri.next();
        }
    }
}



void
starling_pos_processor_base::
pileup_pos_reads(const pos_t pos)
{
    const unsigned sampleCount(getSampleCount());
    for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
    {
        read_segment_iter ri(sample(sampleIndex).readBuffer.get_pos_read_segment_iter(pos));
        read_segment_iter::ret_val r;
        while (true)
        {
            r=ri.get_ptr();
            if (nullptr==r.first) break;
            const read_segment& rseg(r.first->get_segment(r.second));
            pileup_read_segment(rseg, sampleIndex);
            ri.next();
        }
    }
}



void
starling_pos_processor_base::
pileup_read_segment(
    const read_segment& rseg,
    const unsigned sampleIndex)
{
    // get the best alignment for the read:
    const alignment* best_al_ptr(&(rseg.getInputAlignment()));
    if (rseg.is_realigned)
    {
        best_al_ptr=&(rseg.realignment);
    }
    else
    {
        // detect whether this read has no alignments with indels we can handle:
        if (! rseg.is_any_nonovermax(_opt.maxIndelSize)) return;
    }

    const alignment& best_al(*best_al_ptr);
    if (best_al.empty())
    {
        if (! rseg.is_realigned)
        {
            if (_opt.verbosity >= LOG_LEVEL::ALLWARN)
            {
                log_os << "WARNING: skipping read_segment with no genomic alignment and contig alignment outside of indel.\n";
                log_os << "\tread_name: " << rseg.key() << "\n";
            }
        }
        else
        {
            // this indicates that 2 or more equally likely alignments
            // of the entire read exist in the local realignment
            log_os << "WARNING: skipping read_segment which has multiple equally likely but incompatible alignments: " << rseg.key() << "\n";
        }
        return;
    }

    // check that read has not been realigned too far to the left:
    //
    // A warning has already been issued for this at the end of the realignment procedure:
    //
    if (rseg.is_realigned && rseg.is_invalid_realignment) return;

    const unsigned read_size(rseg.read_size());
    const bam_seq bseq(rseg.get_bam_read());
    const uint8_t* qual(rseg.qual());

    // these values are only required for special RNA scoring features:
    const unsigned full_read_size(rseg.full_read_size());
    const unsigned full_read_offset(rseg.full_read_offset());

    static const uint8_t min_adjust_mapq(5);
    const uint8_t mapq(rseg.map_qual());
    const uint8_t adjustedMapq(std::max(min_adjust_mapq,mapq));
    const bool is_mapq_adjust(_opt.isBasecallQualAdjustedForMapq && (adjustedMapq<=80));
    // test read against max indel size (this is a backup, should have been taken care of upstream):
    const unsigned read_ref_mapped_size(apath_ref_length(best_al.path));
    if (read_ref_mapped_size > (read_size+get_largest_total_indel_ref_span_per_read()))
    {
        //brc.large_ref_deletion++;
        return;
    }

    // exact begin and end report range filters:
    {
        if (best_al.pos>=_reportRange.end_pos()) return;
        const pos_t al_end_pos(best_al.pos+static_cast<pos_t>(read_ref_mapped_size));
        if (al_end_pos <= _reportRange.begin_pos()) return;
    }

    /// Find ambiguous sections of the read to trim off
    const unsigned readAmbiguousEndLength(getReadAmbiguousEndLength(bseq, best_al.is_fwd_strand));
    assert(read_size>=readAmbiguousEndLength);

    // read_begin,read_end define a zero-indexed, half-open range in read-coordinates: [read_begin, read_end)
    unsigned read_begin(0);
    unsigned read_end(read_size);

    if (readAmbiguousEndLength > 0)
    {
        if (best_al.is_fwd_strand)
        {
            read_end -= readAmbiguousEndLength;
        }
        else
        {
            read_begin += readAmbiguousEndLength;
        }
    }

    // don't add positions close to the read edge (this is used for error stats estimation)
    if (_opt.minDistanceFromReadEdge > 0)
    {
        read_begin += _opt.minDistanceFromReadEdge;
        if (_opt.minDistanceFromReadEdge <= read_end)
        {
            read_end -= _opt.minDistanceFromReadEdge;
        }
        else
        {
            read_end = 0;
        }

        if (read_end <= read_begin) return;
    }

#ifdef DEBUG_PPOS
    log_os << "best_al: " << best_al;
    log_os << "read_size,read_rms: " << read_size << " " << read_ref_mapped_size << "\n";
    log_os << "read_begin,read_end: " << read_begin << " " << read_end << "\n";
#endif

    // find mismatch filtration:

    // value used for is_neighbor_mismatch in case it is not measured:
    bool is_neighbor_mismatch(false);

    const bool is_submapped(! rseg.is_tier1or2_mapping());
    const bool is_tier1(rseg.is_tier1_mapping());

    // precompute mismatch density info for this read:
    if ((! is_submapped) && _opt.isMismatchDensityFilter())
    {
        const rc_segment_bam_seq ref_bseq(_ref);
        create_mismatch_filter_map(_opt,best_al,ref_bseq,bseq,read_begin,read_end, _candidateSnvBuffer, _rmi);
        if (_opt.useTier2Evidence)
        {
            const int max_pass(_opt.tier2.mismatchDensityFilterMaxMismatchCount);
            for (unsigned i(0); i<read_size; ++i)
            {
                _rmi[i].tier2_mismatch_filter_map = (max_pass < _rmi[i].mismatch_count);
            }
        }
    }

    // info used to determine filter certain deletions, as follows:
    // any edge deletion is thrown away, unless it is next to an exon
    const std::pair<bool,bool> edge_pins(rseg.get_segment_edge_pin());
    const std::pair<unsigned,unsigned> edge_ends(ALIGNPATH::get_match_edge_segments(best_al.path));

    // alignment walkthough:
    pos_t ref_head_pos(best_al.pos);
    unsigned read_head_pos(0);

    // used to phase from the pileup:
    using namespace ALIGNPATH;
    const unsigned as(best_al.path.size());
    for (unsigned i(0); i<as; ++i)
    {
        const path_segment& ps(best_al.path[i]);

#ifdef DEBUG_PPOS
        log_os << "seg,ref,read: " << i << " " << ref_head_pos << " " << read_head_pos << "\n";
#endif

        if (is_segment_align_match(ps.type))
        {
            for (unsigned j(0); j<ps.length; ++j)
            {
                const unsigned read_pos(read_head_pos+j);

                if ((read_pos < read_begin) || (read_pos >= read_end)) continue; // allow for read end trimming
                const pos_t ref_pos(ref_head_pos+static_cast<pos_t>(j));
                // skip position outside of report range:
                if (! is_pos_reportable(ref_pos)) continue;

#ifdef DEBUG_PPOS
                log_os << "j,ref,read: " << j << " " << ref_pos <<  " " << read_pos << "\n";
#endif

                _stagemanPtr->validate_new_pos_value(ref_pos,STAGE::get_pileup_stage_no(_opt));

                const uint8_t call_code(bseq.get_code(read_pos));
                const uint8_t call_id(bam_seq_code_to_id(call_code));

                uint8_t qscore(qual[read_pos]);
                if (is_mapq_adjust)
                {
                    qscore = qphred_to_mapped_qphred(qscore,adjustedMapq);
                }

                const unsigned end_trimmed_read_len(read_end-read_begin);
                unsigned align_strand_read_pos(read_pos);
                if (! best_al.is_fwd_strand)
                {
                    align_strand_read_pos=read_size-(read_pos+1);
                }

                bool current_call_filter( true );
                bool is_tier_specific_filter( false );
                if (! is_submapped)
                {
                    bool is_call_filter((call_code == BAM_BASE::ANY) ||
                                        (qscore < _opt.minBasecallErrorPhredProb));

                    bool is_tier2_call_filter(is_call_filter);
                    if ((! is_call_filter) && _opt.isMismatchDensityFilter())
                    {
                        is_call_filter = _rmi[read_pos].mismatch_filter_map;
                        if (_opt.useTier2Evidence)
                        {
                            is_tier2_call_filter = _rmi[read_pos].tier2_mismatch_filter_map;
                        }
                        else
                        {
                            is_tier2_call_filter = is_call_filter;
                        }
                    }
                    current_call_filter = ( is_tier1 ? is_call_filter : is_tier2_call_filter );
                    is_tier_specific_filter = ( is_tier1 && is_call_filter && (! is_tier2_call_filter) );

                    if (_opt.isMismatchDensityFilter())
                    {
                        is_neighbor_mismatch=(_rmi[read_pos].mismatch_count_ns>0);
                    }
                }

                // always update MAPQ (even when we don't want EVS metrics)
                sample(sampleIndex).basecallBuffer.insert_mapq_count(ref_pos,mapq);

                // update extended feature metrics (including submapped reads):
                if (_opt.is_compute_germline_scoring_metrics())
                {
                    const unsigned full_read_pos(read_pos+full_read_offset);
                    /// zero-indexed distance from the edge of the full read
                    const unsigned distanceFromReadEdge(std::min(full_read_pos, full_read_size-(full_read_pos+1)));

                    updateGermlineScoringMetrics(ref_pos, sampleIndex, call_id, qscore, mapq, align_strand_read_pos,
                                                 distanceFromReadEdge, is_submapped);
                }
                else if (_opt.is_compute_somatic_scoring_metrics)
                {
                    updateSomaticScoringMetrics(ref_pos, sampleIndex, is_tier1, call_id, current_call_filter, read_pos,
                                                read_size);
                }

                if (is_submapped)
                {
                    insert_pos_submap_count(ref_pos,sampleIndex);
                    continue;
                }

                /// include only data meeting mapping criteria after this point:

                try
                {
                    const base_call bc = base_call(call_id,qscore,best_al.is_fwd_strand,
                                                   align_strand_read_pos,end_trimmed_read_len,
                                                   current_call_filter,is_neighbor_mismatch,is_tier_specific_filter);

                    insert_pos_basecall(ref_pos, sampleIndex, is_tier1, bc);
                }
                catch (...)
                {
                    log_os << "Exception caught in starling_pos_processor_base.insert_pos_basecall() "
                           << "while processing read_position: " << (read_pos+1) << "\n";
                    throw;
                }
            }

        }
        else if (ps.type==DELETE)
        {
            const bool is_edge_deletion((i<edge_ends.first)  || (i>edge_ends.second));
            const bool is_pinned_deletion(((i<edge_ends.first) && (edge_pins.first)) ||
                                          ((i>edge_ends.second) && (edge_pins.second)));
            if ((! is_edge_deletion) || is_pinned_deletion)
            {
                for (unsigned j(0); j<ps.length; ++j)
                {
                    const pos_t ref_pos(ref_head_pos+static_cast<pos_t>(j));

                    // skip position outside of report range:
                    if (! is_pos_reportable(ref_pos)) continue;
                    _stagemanPtr->validate_new_pos_value(ref_pos,STAGE::get_pileup_stage_no(_opt));

                    if (is_submapped)
                    {
                        insert_pos_submap_count(ref_pos,sampleIndex);
                    }
                    else
                    {
                        insert_pos_spandel_count(ref_pos,sampleIndex);
                    }
                }
            }
        }

        if (is_segment_type_read_length(ps.type)) read_head_pos += ps.length;
        if (is_segment_type_ref_length(ps.type)) ref_head_pos += ps.length;
    }

    //    return READ_FATE::USED;
}



void
starling_pos_processor_base::
process_pos_variants(
    const pos_t pos,
    const bool isPosPrecedingReportableRange)
{
    if (not isPosPrecedingReportableRange)
    {
        process_pos_stats(pos);

        const unsigned sampleCount(getSampleCount());
        for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
        {
            process_pos_sample_stats(pos, sampleIndex);
        }
    }
    process_pos_variants_impl(pos, isPosPrecedingReportableRange);
}



void
starling_pos_processor_base::
process_pos_stats(
    const pos_t pos)
{
    {
        // update indel stats
        auto it(getIndelBuffer().positionIterator(pos));
        const auto it_end(getIndelBuffer().positionIterator(pos + 1));

        for (; it != it_end; ++it)
        {
            const IndelKey& indelKey(it->first);
            const IndelData& indelData(getIndelData(it));

            if (indelKey.isMismatch()) continue;

            const bool isCandidate(getIndelBuffer().isCandidateIndel(indelKey, indelData));
            _statsManager.addCallRegionIndel(isCandidate);
        }
    }
}



void
starling_pos_processor_base::
process_pos_sample_stats(
    const pos_t pos,
    const unsigned sample_no)
{
    sample_info& sif(sample(sample_no));

    const snp_pos_info& pi(sif.basecallBuffer.get_pos(pos));

    static const bool is_include_tier2(false);
    _pileupCleaner.CleanPileupFilter(pi,is_include_tier2,sif.cleanedPileup);

    const unsigned n_spandel(pi.spanningDeletionReadCount);
    const unsigned n_submapped(pi.submappedReadCount);

    if (pi.get_ref_base() != 'N')
    {
        sif.localRegionStatsCollection.insert(pos, sif.cleanedPileup.usedBasecallCount(), sif.cleanedPileup.unusedBasecallCount(),n_spandel,n_submapped);
    }
    else
    {
        sif.localRegionStatsCollection.insert_null(pos);
    }
}

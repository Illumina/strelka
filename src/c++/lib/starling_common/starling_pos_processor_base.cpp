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
///

#ifdef _MSC_VER
#pragma warning(disable:4355)
#endif

#include "starling_pos_processor_indel_util.hh"
#include "starling_read_align.hh"
#include "starling_read_util.hh"

#include "blt_common/position_snp_call_pprob_digt.hh"
#include "blt_common/position_snp_call_pprob_nploid.hh"
#include "blt_util/depth_buffer_util.hh"
#include "blt_util/io_util.hh"
#include "blt_util/log.hh"
#include "htsapi/bam_seq_read_util.hh"
#include "starling_common/starling_indel_report_info.hh"

#include <iomanip>


//#define DEBUG_PPOS

// largest_read_size grows dynamically with observed read size, this
// is used to initialize the setting prior to observing any reads:
//
// initial setting is large to help consistently deal with grouperisms:
//
const unsigned STARLING_INIT_LARGEST_READ_SIZE(250);
const double STARLING_LARGEST_READ_SIZE_PAD(1.25);

// largest indel_size grows dynamically with observed indel size until
// hitting max_indel_size. Initialized to the follow value prior to
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



#if 0
static
void
write_snp_prefix_info(const char* label,
                      const pos_t output_pos,
                      const char ref,
                      const unsigned n_used_calls,
                      const unsigned n_unused_calls,
                      std::ostream& os)
{

    os << label
       << " pos: " << output_pos
       << " bcalls_used: " << n_used_calls
       << " bcalls_filt: " << n_unused_calls
       << " ref: " << ref;
}
#endif



static
unsigned
get_read_buffer_size(const unsigned largest_read_size,
                     const unsigned largest_total_indel_span_per_read)
{
    return (largest_read_size+largest_total_indel_span_per_read);
}



// public companion functions:
//
static
int
get_influence_zone_size(const unsigned largest_read_size,
                        const unsigned largest_total_indel_span_per_read)
{
    static const unsigned min_influence_zone_read_size(512);
    const unsigned influence_read_size(std::max(min_influence_zone_read_size,
                                                largest_read_size));
    return static_cast<int>(get_read_buffer_size(influence_read_size,largest_total_indel_span_per_read))-1;
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
    sdata.add_stage(READ_BUFFER,HEAD,get_read_buffer_size(largest_read_size,largest_total_indel_ref_span_per_read));

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
starling_pos_processor_base(const starling_base_options& opt,
                            const starling_base_deriv_options& dopt,
                            const reference_contig_segment& ref,
                            const starling_streams_base& streams,
                            const unsigned n_samples)
    : base_t()
    , _opt(opt)
    , _dopt(dopt)
    , _ref(ref)
    , _streams(streams)
    , _rmi(STARLING_INIT_LARGEST_READ_SIZE)
      //, _largest_indel_size(std::min(opt.max_indel_size,STARLING_INIT_LARGEST_INDEL_SIZE)) -- tmp change for GRUOPER handling
    , _largest_indel_ref_span(opt.max_indel_size)
    , _largest_total_indel_ref_span_per_read(_largest_indel_ref_span)
    , _stageman(STAGE::get_stage_data(STARLING_INIT_LARGEST_READ_SIZE, get_largest_total_indel_ref_span_per_read(), _opt, _dopt),dopt.report_range,*this)
    , _chrom_name(_opt.bam_seq_name)
    , _n_samples(n_samples)
    , _pileupCleaner(opt)
    , _indelBuffer(opt,dopt,ref)
    , _active_region_detector(ref, _indelBuffer)
{
    assert((_n_samples != 0) && (_n_samples <= MAX_SAMPLE));

    const unsigned report_size(_dopt.report_range.size());
    const unsigned knownref_report_size(get_ref_seq_known_size(_ref,_dopt.report_range));
    for (unsigned i(0); i<_n_samples; ++i)
    {
        _sample[i].reset(new sample_info(_opt, ref, report_size,knownref_report_size,&_ric));
    }

#ifdef HAVE_FISHER_EXACT_TEST
    if (_opt.is_adis_table)
    {
        _ws.reset(get_exact_test_ws());
    }
#endif

    if (_opt.is_bsnp_nploid)
    {
        _ninfo.reset(new nploid_info(_opt.bsnp_nploid_ploidy));
    }

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

    // define an expanded indel influence zone around the report range:
    //
    // note that we don't know the max indel ref span per read at this point, so a fudge factor is
    // added here:
    const int bshift(get_influence_zone_size(get_largest_read_size(),
                                             _opt.max_indel_size*2));
    pos_range& rir( _report_influence_range);
    rir = _dopt.report_range_limit;
    if (rir.is_begin_pos)
    {
        rir.begin_pos -= bshift;
    }
    if (rir.is_end_pos)
    {
        rir.end_pos += bshift;
    }
}



void
starling_pos_processor_base::
update_stageman()
{
    _stageman.revise_stage_data(
        STAGE::get_stage_data(get_largest_read_size(),
                              get_largest_total_indel_ref_span_per_read(),
                              _opt,
                              _dopt));
}


bool
starling_pos_processor_base::
update_largest_read_size(const unsigned rs)
{
    if (rs>STARLING_MAX_READ_SIZE) return false;

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
    assert(is<=_opt.max_indel_size);
    _largest_indel_ref_span=std::min(is,_opt.max_indel_size);
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
    _stageman.reset();

    pos_range output_report_range(_dopt.report_range);

    if ((! output_report_range.is_begin_pos) &&
        _stageman.is_first_pos_set())
    {
        output_report_range.set_begin_pos(_stageman.min_pos());
    }

    if ((! output_report_range.is_end_pos) &&
        _stageman.is_first_pos_set())
    {
        output_report_range.set_end_pos(_stageman.max_pos()+1);
    }

    write_counts(output_report_range);
}



bool
starling_pos_processor_base::
insert_indel(
    const IndelObservation& obs,
    const unsigned sampleId)
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
        _stageman.validate_new_pos_value(obs.key.pos,STAGE::READ_BUFFER);

        // dynamically scale maximum indel size:
        const unsigned len(std::min(static_cast<unsigned>((obs.key.delete_length())),_opt.max_indel_size));
        update_largest_indel_ref_span(len);

        bool is_novel(getIndelBuffer().addIndelObservation(sampleId, obs));
        if (obs.data.is_forced_output) _is_skip_process_pos=false;

        return is_novel;
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
    _stageman.validate_new_pos_value(pos,STAGE::READ_BUFFER);
    _forced_output_pos.insert(pos);
    _is_skip_process_pos=false;
}



bool
starling_pos_processor_base::
insert_ploidy_region(
    const known_pos_range2& range,
    const unsigned ploidy)
{
    assert(ploidy==0 || ploidy==1);
    _stageman.validate_new_pos_value(range.begin_pos(),STAGE::READ_BUFFER);
    return _ploidy_regions.addRegion(range,ploidy);
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


// TODO use boost::optional here:
//
boost::optional<align_id_t>
starling_pos_processor_base::
insert_read(
    const bam_record& br,
    const alignment& al,
    const char* chrom_name,
    const MAPLEVEL::index_t maplev,
    const unsigned sample_no)
{
    boost::optional<align_id_t> retval;

    if (0 != strcmp(_chrom_name.c_str(),chrom_name))
    {
        log_os << "ERROR: starling_pos_processor_base.insert_read(): read has unexpected sequence name: '" << chrom_name << "' expecting: '" << _chrom_name << "'\n"
               << "\tread_key: " << read_key(br) << "\n";
        exit(EXIT_FAILURE);
    }

    // Check at al.pos rather than buffer pos because we need to
    // insure that all candidate indels get entered by the end of HEAD
    // stage. If the read aligns back past the end of head stage it is
    // ok so long as we know it will not generate any candidate indels
    // in this region:
    //
    if (! _stageman.is_new_pos_value_valid(al.pos,STAGE::HEAD))
    {
        log_os << "WARNING: skipping alignment for read: " << read_key(br)
               << " which falls outside of the read buffer\n";
        return retval;
    }

    starling_read_buffer& rbuff(sample(sample_no).read_buff);

    // check whether the read buffer has reached max capacity
    if (_opt.isMaxBufferedReads())
    {
        if (rbuff.size() >= _opt.maxBufferedReads)
        {
            return retval;
        }
    }

    // assume that pos_procesor, as a container, is no longer empty...
    _is_skip_process_pos=false;

    // update read_size:
    {
        const unsigned rs(static_cast<unsigned>(STARLING_LARGEST_READ_SIZE_PAD*br.read_size()));
        if (! update_largest_read_size(rs))
        {
            std::ostringstream oss;
            oss << "ERROR: Input read size: " << br.read_size() << " exceeds maximum.";
            throw blt_exception(oss.str().c_str());
        }
    }

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
        load_read_in_depth_buffer(sread_ptr->get_full_segment(),sample_no);

        // update other data for only the first read segment
        const seg_id_t seg_id(sread_ptr->is_segmented() ? 1 : 0 );
        init_read_segment(sread_ptr->get_segment(seg_id),sample_no);
    }

    return retval;
}



// // this function is last chance to check-for/warn-about/clean-out any
// // alignments that won't survive alignment and calling:
// //
// // note that cleared alignments still potentially have ids sitting in
// // the indel buffer, and this would be hard to clean out (except by
// // brute force)
// //
// void
// starling_pos_processor_base::
// clean_pos(const pos_t pos) {
//     std::vector<align_id_t> dead_meat;

//     starling_read_iter ri(_read_buff.get_pos_read_iter(pos));
//     starling_read* srp;
//     while(NULL!=(srp=ri.get_ptr())){
//         if(not srp->is_full_record()) {
//             log_os << "WARNING: incomplete read record must be removed from pipeline: " << srp->key() << "\n";
//             dead_meat.push_back(srp->id);
//         }
//         ri.next();
//     }

//     const unsigned ds(dead_meat.size());
//     for(unsigned i(0);i<ds;++i) _read_buff.clear_read(dead_meat[i]);
// }



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
    const alignment& al(rseg.genome_align());
    if (al.empty()) return;

    const MAPLEVEL::index_t maplev(rseg.genome_align_maplev());
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



// only acts on genomic mapped reads:
void
starling_pos_processor_base::
init_read_segment(
    const read_segment& rseg,
    const unsigned sample_no)
{
    const alignment& al(rseg.genome_align());
    if (al.empty()) return;

    const MAPLEVEL::index_t maplev(rseg.genome_align_maplev());
    const INDEL_ALIGN_TYPE::index_t iat(translate_maplev_to_indel_type(maplev));

    const bam_seq bseq(rseg.get_bam_read());
    try
    {
        const unsigned total_indel_ref_span_per_read =
            add_alignment_indels_to_sppr(_opt.max_indel_size,_ref,
                                         al,bseq,*this,iat,rseg.id(),sample_no,rseg.get_segment_edge_pin(), rseg.map_qual() == 0);

        update_largest_total_indel_ref_span_per_read(total_indel_ref_span_per_read);

    }
    catch (...)
    {
        log_os << "\nException caught in add_alignment_indels_to_sppr() while processing record: " << rseg << "\n";
        throw;
    }
}



// For all read segments buffered at the current position:
// 1) process genomic alignment of read segment for indels
// 2) add genomic alignment of read segment to estdepth
//
void
starling_pos_processor_base::
init_read_segment_pos(const pos_t pos)
{
    for (unsigned s(0); s<_n_samples; ++s)
    {
        read_segment_iter ri(sample(s).read_buff.get_pos_read_segment_iter(pos));
        for (read_segment_iter::ret_val r; true; ri.next())
        {
            r=ri.get_ptr();
            if (NULL==r.first) break;
            // full_segments of unspliced reads and the initial
            // segment of spliced reads are initialized outside of the
            // process_pos framework, so this routine only initializes
            // the second segment or higher:
            //
            if (r.second<2) continue;
            init_read_segment(r.first->get_segment(r.second),s);

        }
    }
}



#if 0
// function has the general purpose of normalizing candidates which
// are very likely to refer to the same event -- in practice at the
// moment this only includes removing breakpoint calls which fall
// entirely within an existing closed insertion
//
void
consolidate_candidate_indel_pos(pos)
{


}
#endif



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



// For all reads buffered at the current position:
// 1) determine the set of candidate indels that the read overlaps
// 2) determine the set of private indels within the read's discovery alignments
// 3) Find the most likely alignment given both sets of indels
// 4) evaluate the probability that the read supports each candidate indel vs. the reference
//
// these occur in subsequent methods:
//
// 5) process most-likely alignment for snp-calling
// 5a) insert most-likely alignment into output read buffer (re-link within same data-structure?)
// 6) remove read from input read buffer

void
starling_pos_processor_base::
align_pos(const pos_t pos)
{
    known_pos_range realign_buffer_range(get_realignment_range(pos, _stageman.get_stage_data()));

    for (unsigned sampleIndex(0); sampleIndex<_n_samples; ++sampleIndex)
    {
        sample_info& sif(sample(sampleIndex));
        read_segment_iter ri(sif.read_buff.get_pos_read_segment_iter(pos));
        for (read_segment_iter::ret_val r; true; ri.next())
        {
            r=ri.get_ptr();
            if (nullptr == r.first) break;
            read_segment& rseg(r.first->get_segment(r.second));
            if (_opt.is_realign_submapped_reads ||
                rseg.is_tier1or2_mapping())
            {
                try
                {
                    realign_and_score_read(_opt,_dopt,sif.sample_opt,_ref,realign_buffer_range,sampleIndex,rseg,
                                           getIndelBuffer());
                }
                catch (...)
                {
                    log_os << "ERROR: Exception caught in align_pos() while realigning segment: "
                           << static_cast<int>(r.second) << " of read: " << (*r.first) << "\n";
                    throw;
                }
                // check that read has not been realigned too far to the left:
                if (rseg.is_realigned)
                {
                    if (! _stageman.is_new_pos_value_valid(rseg.realignment.pos,STAGE::POST_ALIGN))
                    {
                        log_os << "WARNING: read realigned outside bounds of realignment stage buffer. Skipping...\n"
                               << "\tread: " << rseg.key() << "\n";
                        rseg.is_invalid_realignment=true;
                    }
                }
            }
        }

#ifdef ESTDEPTH_DUMP
        if (sif.estdepth_buff.val(pos)>0)
        {
            log_os << "ESTDEPTH: pos,val: " << pos << " " << sif.estdepth_buff.val(pos) << "\n";
        }
#endif
    }
}



void
starling_pos_processor_base::
set_head_pos(const pos_t pos)
{
    _stageman.validate_new_pos_value(pos,STAGE::READ_BUFFER);
    _stageman.handle_new_pos_value(pos);
}



void
starling_pos_processor_base::
process_pos(const int stage_no,
            const pos_t pos)
{
#if 0
    log_os << "pos,stage_no: " << pos << " " << stage_no << "\n";
#endif

    if (empty()) return;

    if        (stage_no==STAGE::HEAD)
    {
        init_read_segment_pos(pos);
        if (_opt.is_short_haplotype_calling_enabled)
        {
            _active_region_detector.updateEndPosition(pos);

        }

        if (_opt.is_write_candidate_indels())
        {
            if (is_pos_reportable(pos))
            {
                write_candidate_indels_pos(pos);
            }
        }
    }
    else if (stage_no==STAGE::READ_BUFFER)
    {
#if 0
        for (unsigned s(0); s<_n_samples; ++s)
        {
            sample_info& sif(sample(s));
            sif.indel_buff.dump_pos(pos,log_os);
            sif.read_buff.dump_pos(pos,log_os);
        }
#endif
        //        consolidate_candidate_indel_pos(pos);
        if (_opt.is_short_haplotype_calling_enabled)
            _active_region_detector.updateStartPosition(pos);

        if (! _opt.is_write_candidate_indels_only)
        {
            //        clean_pos(pos);
            align_pos(pos);
            pileup_pos_reads(pos);
            // if(_opt.is_realigned_read_file) {
            //     rebuffer_pos_reads(pos);
            // }

            write_reads(pos);
        }

    }
    else if (stage_no==STAGE::POST_ALIGN)
    {
        if (! _opt.is_write_candidate_indels_only)
        {
            if (is_pos_reportable(pos))
            {
                process_pos_variants(pos);
            }
        }

        // everything else:
        post_align_clear_pos(pos);
    }
    else if (stage_no==STAGE::CLEAR_SITE_ANNOTATION)
    {
        _forced_output_pos.erase(pos);
        _ploidy_regions.removeToPos(pos);
        clear_pos_annotation(pos);
    }
    else if (stage_no==STAGE::CLEAR_READ_BUFFER)
    {
        for (unsigned s(0); s < _n_samples; ++s)
        {
            sample(s).read_buff.clear_to_pos(pos);
        }
    }
    else if (stage_no==STAGE::CLEAR_INDEL_BUFFER)
    {
        getIndelBuffer().clearPosition(pos);
    }
    else if (stage_no==STAGE::POST_CALL)
    {
        for (unsigned s(0); s<_n_samples; ++s)
        {
            sample_info& sif(sample(s));
            sif.estdepth_buff.clear_pos(pos);
            sif.estdepth_buff_tier2.clear_pos(pos);

            // if we are doing short-range phasing, suspend pileup buffer clear
            // while phasing block is being built:
            if (! is_save_pileup_buffer())
            {
                sif.bc_buff.clear_to_pos(pos);
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

    sample(sample_no).bc_buff.insert_pos_submap_count(pos);
}



void
starling_pos_processor_base::
insert_pos_spandel_count(const pos_t pos,
                         const unsigned sample_no)
{
    // assume pos has been pre-checked:

    sample(sample_no).bc_buff.insert_pos_spandel_count(pos);
}



void
starling_pos_processor_base::
update_ranksum_and_mapq_count(
    const pos_t pos,
    const unsigned sample_no,
    const uint8_t call_id,
    const uint8_t qscore,
    const uint8_t mapq,
    const unsigned cycle,
    const bool is_submapped)
{
    // assume pos is already valid:

    auto& bcbuff(sample(sample_no).bc_buff);
    bcbuff.insert_mapq_count(pos,mapq);
    bcbuff.update_ranksums(_ref.get_base(pos),pos,call_id,qscore,mapq,cycle,is_submapped);
}


void
starling_pos_processor_base::
update_somatic_features(
    const pos_t pos,
    const unsigned sample_no,
    const bool is_tier1,
    const uint8_t call_id,
    const bool is_call_filter,
    const uint8_t mapq,
    const uint16_t readPos,
    const uint16_t readLength)
{
    // assume pos is already valid:

    auto& bcbuff(sample(sample_no).bc_buff);
    bcbuff.insert_mapq_count(pos,mapq);

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

    sample(sample_no).bc_buff.insert_pos_basecall(pos,is_tier1,bc);
}



void
starling_pos_processor_base::
insert_hap_cand(const pos_t pos,
                const unsigned sample_no,
                const bool is_tier1,
                const bam_seq_base& read_seq,
                const uint8_t* qual,
                const unsigned offset)
{
    // assume pos is already valid:

    sample(sample_no).bc_buff.insert_hap_cand(pos,is_tier1,read_seq,qual,offset);
}



void
starling_pos_processor_base::
write_candidate_indels_pos(
    const pos_t pos)
{
    const pos_t output_pos(pos+1);
    std::ostream& bos(*_streams.candidate_indel_osptr());

    auto indelIter(getIndelBuffer().positionIterator(pos));
    const auto indelIterEnd(getIndelBuffer().positionIterator(pos + 1));

    for (; indelIter!=indelIterEnd; ++indelIter)
    {
        const IndelKey& indelKey(indelIter->first);
        const IndelData& indelData(getIndelData(indelIter));
        if (!getIndelBuffer().isCandidateIndel(indelKey, indelData)) continue;
        bos << _chrom_name << "\t"
            << output_pos << "\t"
            << INDEL::get_index_label(indelKey.type) << "\t";
        if (indelKey.is_breakpoint())
        {
            bos << 0 << "\t" << indelData.getBreakpointInsertSeq();
        }
        else
        {
            bos << indelKey.deletionLength << "\t"
                << indelKey.insert_seq();
        }
        bos << "\n";
    }
}



#if 1
static
pos_t
get_new_read_pos(const read_segment& rseg)
{
    // get the best alignment for the read:
    const alignment* best_al_ptr(&(rseg.genome_align()));
    if (rseg.is_realigned) best_al_ptr=&(rseg.realignment);

    if (best_al_ptr->empty()) return rseg.buffer_pos;     // a grouper contig read which was not realigned...
    else                     return best_al_ptr->pos;
}



// adjust read buffer position so that reads are buffered in sorted
// order after realignment:
//
void
starling_pos_processor_base::
rebuffer_pos_reads(const pos_t pos)
{
    // need to queue up read changes and run at the end so that we
    // don't invalidate read buffer iterators
    //
    typedef std::pair<std::pair<align_id_t,seg_id_t>,pos_t> read_pos_t;

    for (unsigned s(0); s<_n_samples; ++s)
    {
        sample_info& sif(sample(s));
        std::vector<read_pos_t> new_read_pos;
        read_segment_iter ri(sif.read_buff.get_pos_read_segment_iter(pos));
        for (read_segment_iter::ret_val r; true; ri.next())
        {
            r=ri.get_ptr();
            if (NULL==r.first) break;
            read_segment& rseg(r.first->get_segment(r.second));

            const pos_t new_pos(get_new_read_pos(rseg));
            if ((new_pos!=pos) &&
                (_stageman.is_new_pos_value_valid(new_pos,STAGE::POST_ALIGN)))
            {
                new_read_pos.push_back(std::make_pair(std::make_pair(rseg.id(),r.second),new_pos));
            }
        }

        const unsigned nr(new_read_pos.size());
        for (unsigned i(0); i<nr; ++i)
        {
            sif.read_buff.rebuffer_read_segment(new_read_pos[i].first.first,
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
    for (unsigned s(0); s<_n_samples; ++s)
    {
        bam_dumper* bamd_ptr(_streams.realign_bam_ptr(s));
        if (NULL == bamd_ptr) continue;
        bam_dumper& bamd(*bamd_ptr);

        read_segment_iter ri(sample(s).read_buff.get_pos_read_segment_iter(pos));
        read_segment_iter::ret_val r;

        while (true)
        {
            r=ri.get_ptr();
            if (nullptr==r.first) break;
            if (r.first->segment_count()==r.second)
            {
                r.first->write_bam(bamd);
            }
            ri.next();
        }
    }
}



// convert reads buffered at position into a position basecall
// "pileup" to allow for downstream depth and snp calculations
//
void
starling_pos_processor_base::
pileup_pos_reads(const pos_t pos)
{
    for (unsigned s(0); s<_n_samples; ++s)
    {
        read_segment_iter ri(sample(s).read_buff.get_pos_read_segment_iter(pos));
        read_segment_iter::ret_val r;
        while (true)
        {
            r=ri.get_ptr();
            if (nullptr==r.first) break;
            const read_segment& rseg(r.first->get_segment(r.second));
            pileup_read_segment(rseg,s);
            ri.next();
        }
    }
}



void
starling_pos_processor_base::
pileup_read_segment(const read_segment& rseg,
                    const unsigned sample_no)
{
    // get the best alignment for the read:
    const alignment* best_al_ptr(&(rseg.genome_align()));
    if (rseg.is_realigned)
    {
        best_al_ptr=&(rseg.realignment);
    }
    else
    {
        // detect whether this read has no alignments with indels we can handle:
        if (! rseg.is_any_nonovermax(_opt.max_indel_size)) return;
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
    // A warning has already been issues for this at the end of realignment:
    //
    if (rseg.is_realigned && rseg.is_invalid_realignment) return;

    const unsigned read_size(rseg.read_size());
    const bam_seq bseq(rseg.get_bam_read());
    const uint8_t* qual(rseg.qual());

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
        const pos_range& rlimit(_dopt.report_range_limit);
        if (rlimit.is_end_pos && (best_al.pos>=rlimit.end_pos)) return;
        if (rlimit.is_begin_pos)
        {
            const pos_t al_end_pos(best_al.pos+static_cast<pos_t>(read_ref_mapped_size));
            if (al_end_pos <= rlimit.begin_pos) return;
        }
    }

    // find trimmed sections (as defined by the CASAVA 1.0 caller)
    unsigned fwd_strand_begin_skip(0);
    unsigned fwd_strand_end_skip(0);
    get_read_fwd_strand_skip(bseq,
                             best_al.is_fwd_strand,
                             fwd_strand_begin_skip,
                             fwd_strand_end_skip);

    assert(read_size>=fwd_strand_end_skip);
    const unsigned read_begin(fwd_strand_begin_skip);
    const unsigned read_end(read_size-fwd_strand_end_skip);

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
    if ((! is_submapped) && _opt.is_max_win_mismatch)
    {
        const rc_segment_bam_seq ref_bseq(_ref);
        create_mismatch_filter_map(_opt,best_al,ref_bseq,bseq,read_begin,read_end,_active_region_detector, _rmi);
        if (_opt.tier2.is_tier2_mismatch_density_filter_count)
        {
            const int max_pass(_opt.tier2.tier2_mismatch_density_filter_count);
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

//                bool isPolySite = _active_region_detector.isPolymorphicSite(ref_pos);

                bool isFirstBaseCallFromMatchSeg(false);
                const bool isFirstFromMatchSeg((j==0) || (read_pos==read_begin) || (! is_pos_reportable(ref_pos-1)));
                if (isFirstFromMatchSeg)
                {
                    if (i==0)
                    {
                        isFirstBaseCallFromMatchSeg=true;
                    }
                    else
                    {
                        const path_segment& last_ps(best_al.path[i-1]);
                        if (! is_segment_align_match(last_ps.type)) isFirstBaseCallFromMatchSeg=true;
                    }

                }

                bool isLastBaseCallFromMatchSeg(false);
                const bool isLastFromMatchSeg(((j+1) == ps.length) || ((read_pos+1) == read_end) || (! is_pos_reportable(ref_pos+1)));
                if (isLastFromMatchSeg)
                {
                    if ((i+1) == as)
                    {
                        isLastBaseCallFromMatchSeg=true;
                    }
                    else
                    {
                        const path_segment& next_ps(best_al.path[i+1]);
                        if (! is_segment_align_match(next_ps.type)) isLastBaseCallFromMatchSeg=true;
                    }
                }

#ifdef DEBUG_PPOS
                log_os << "j,ref,read: " << j << " " << ref_pos <<  " " << read_pos << "\n";
#endif

                _stageman.validate_new_pos_value(ref_pos,STAGE::get_pileup_stage_no(_opt));

                const uint8_t call_code(bseq.get_code(read_pos));
                const uint8_t call_id(bam_seq_code_to_id(call_code));

                uint8_t qscore(qual[read_pos]);
                if (is_mapq_adjust)
                {
                    qscore = qphred_to_mapped_qphred(qscore,adjustedMapq);
                }

                unsigned align_strand_read_pos(read_pos);
                unsigned end_trimmed_read_len(read_end);
                if (! best_al.is_fwd_strand)
                {
                    align_strand_read_pos=read_size-(read_pos+1);
                    end_trimmed_read_len=read_size-fwd_strand_begin_skip;
                }

                bool current_call_filter( true );
                bool is_tier_specific_filter( false );
                if (! is_submapped)
                {
                    bool is_call_filter((call_code == BAM_BASE::ANY) ||
                                        (qscore < _opt.min_qscore));

                    assert(! _opt.is_min_win_qscore);

                    bool is_tier2_call_filter(is_call_filter);
                    if (! is_call_filter && _opt.is_max_win_mismatch)
                    {
                        is_call_filter = _rmi[read_pos].mismatch_filter_map; // && !isPolySite;
                        if (! _opt.tier2.is_tier2_no_mismatch_density_filter)
                        {
                            if (_opt.tier2.is_tier2_mismatch_density_filter_count)
                            {
                                is_tier2_call_filter = _rmi[read_pos].tier2_mismatch_filter_map; // && !isPolySite;
                            }
                            else
                            {
                                is_tier2_call_filter = is_call_filter;
                            }
                        }
                    }
                    current_call_filter = ( is_tier1 ? is_call_filter : is_tier2_call_filter );
                    is_tier_specific_filter = ( is_tier1 && is_call_filter && (! is_tier2_call_filter) );

                    if (_opt.is_max_win_mismatch)
                    {
                        is_neighbor_mismatch=(_rmi[read_pos].mismatch_count_ns>0);
                    }
                }

                // update extended feature metrics (including submapped reads):
                if (_opt.is_compute_germline_scoring_metrics())
                {
                    update_ranksum_and_mapq_count(ref_pos,sample_no,call_id,qscore,mapq,align_strand_read_pos,is_submapped);
                }
                else if (_opt.is_compute_somatic_scoring_metrics)
                {
                    update_somatic_features(ref_pos,sample_no,is_tier1,call_id,current_call_filter,mapq,read_pos,read_size);
                }

                if (is_submapped)
                {
                    insert_pos_submap_count(ref_pos,sample_no);
                    continue;
                }

                /// include only data meeting mapping criteria after this point:

                if (_opt.is_compute_hapscore)
                {
                    insert_hap_cand(ref_pos,sample_no,is_tier1,
                                    bseq,qual,read_pos);
                }

                try
                {
                    const base_call bc = base_call(call_id,qscore,best_al.is_fwd_strand,
                                                   align_strand_read_pos,end_trimmed_read_len,
                                                   current_call_filter,is_neighbor_mismatch,is_tier_specific_filter,
                                                   isFirstBaseCallFromMatchSeg,
                                                   isLastBaseCallFromMatchSeg);

                    insert_pos_basecall(ref_pos,
                                        sample_no,
                                        is_tier1,
                                        bc);
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
                    _stageman.validate_new_pos_value(ref_pos,STAGE::get_pileup_stage_no(_opt));

                    if (is_submapped)
                    {
                        insert_pos_submap_count(ref_pos,sample_no);
                    }
                    else
                    {
                        insert_pos_spandel_count(ref_pos,sample_no);
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
process_pos_variants(const pos_t pos)
{
    for (unsigned sampleId(0); sampleId<_n_samples; ++sampleId)
    {
        process_pos_site_stats(pos,sampleId);
    }
    process_pos_variants_impl(pos);
}



void
starling_pos_processor_base::
process_pos_site_stats(
    const pos_t pos,
    const unsigned sample_no)
{
    sample_info& sif(sample(sample_no));

    const snp_pos_info& pi(sif.bc_buff.get_pos(pos));

    static const bool is_include_tier2(false);
    _pileupCleaner.CleanPileupFilter(pi,is_include_tier2,sif.cpi);

    const unsigned n_spandel(pi.n_spandel);
    const unsigned n_submapped(pi.n_submapped);

    sif.ss.update(sif.cpi.n_calls());
    sif.used_ss.update(sif.cpi.n_used_calls());
    if (pi.get_ref_base() != 'N')
    {
        sif.ssn.update(sif.cpi.n_calls());
        sif.used_ssn.update(sif.cpi.n_used_calls());
        sif.wav.insert(pos,sif.cpi.n_used_calls(),sif.cpi.n_unused_calls(),n_spandel,n_submapped);
    }
    else
    {
        sif.wav.insert_null(pos);
    }
}



const diploid_genotype&
starling_pos_processor_base::
get_empty_dgt(const char ref) const
{
    if (ref!='N')
    {
        return static_cast<const diploid_genotype&>(*_empty_dgt[base_to_id(ref)]);
    }
    else
    {
        static const diploid_genotype n_dgt;
        return static_cast<const diploid_genotype&>(n_dgt);
    }
}

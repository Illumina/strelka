// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Starka
// Copyright (c) 2009-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

///
/// \author Chris Saunders
///

#ifdef _MSC_VER
#pragma warning(disable:4355)
#endif

#include "depth_buffer_util.hh"
#include "starling_pos_processor_indel_util.hh"
#include "starling_read_align.hh"
#include "starling_read_util.hh"

#include "blt_common/adjust_joint_eprob.hh"
#include "blt_common/position_nonref_test.hh"
#include "blt_common/position_nonref_2allele_test.hh"
#include "blt_common/position_snp_call_lrt.hh"
#include "blt_common/position_snp_call_pprob_digt.hh"
#include "blt_common/position_snp_call_pprob_monogt.hh"
#include "blt_common/position_snp_call_pprob_nploid.hh"
#include "blt_common/position_strand_coverage_anomaly.hh"
#include "blt_common/position_strand_distro_anomaly.hh"
#include "blt_common/position_strand_distro_anomaly_lrt.hh"
#include "blt_common/ref_context.hh"
#include "blt_util/align_path.hh"
#include "blt_util/blt_exception.hh"
#include "blt_util/log.hh"
#include "blt_util/read_util.hh"
#include "starling_common/starling_indel_call_pprob_digt.hh"
#include "starling_common/starling_indel_error_prob.hh"
#include "starling_common/starling_indel_report_info.hh"
#include "starling_common/starling_pos_processor_base.hh"
#include "starling_common/gvcf_aggregator.hh"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>


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
report_counts(const snp_pos_info& pi,
              const unsigned n_unused_calls,
              const pos_t output_pos,
              std::ostream& os)
{
    unsigned base_count[N_BASE];

    for (unsigned i(0); i<N_BASE; ++i) base_count[i] = 0;

    for (const auto& bc : pi.calls)
    {
        assert(bc.base_id!=BASE_ID::ANY);
        base_count[bc.base_id]++;
    }

    os << output_pos << '\t';
    for (unsigned i(0); i<N_BASE; ++i)
    {
        os << base_count[i] << '\t';
    }
    os << n_unused_calls << '\n';
}



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
void
write_snp_prefix_info_file(const std::string& seq_name,
                           const pos_t output_pos,
                           const char ref,
                           const unsigned n_used_calls,
                           const unsigned n_unused_calls,
                           std::ostream& os)
{
    os << seq_name << "\t"
       << output_pos << "\t"
       << n_used_calls << "\t"
       << n_unused_calls << "\t"
       << ref;
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
void
write_bsnp_diploid_allele(const blt_options& client_opt,
                          const blt_streams& client_io,
                          const std::string& seq_name,
                          const pos_t output_pos,
                          const char ref,
                          const unsigned n_used_calls,
                          const unsigned n_unused_calls,
                          const snp_pos_info& good_pi,
                          const diploid_genotype& dgt,
                          const unsigned hpol = 0)
{
    std::ostream& os(*client_io.bsnp_diploid_allele_osptr());

    write_snp_prefix_info_file(seq_name,output_pos,ref,n_used_calls,n_unused_calls,os);
    os << "\t";
    write_diploid_genotype_allele(client_opt,good_pi,dgt,os,hpol);
    os << "\n";
}



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
get_pileup_stage_no(const starling_options& opt)
{
    return (opt.is_htype_calling ?
            static_cast<int>(POST_REGION) :
            static_cast<int>(POST_ALIGN));
}

// stage into which pileup entries must fit:
static
int
get_last_static_stage_no(const starling_options& opt)
{
    return (opt.is_htype_calling ?
            static_cast<int>(POST_CALL) :
            static_cast<int>(POST_CALL));
}

// given largest read, and indel ref span per read, get stage
// lengths:
//
static
stage_data
get_stage_data(const unsigned largest_read_size,
               const unsigned largest_total_indel_ref_span_per_read,
               const starling_options& opt,
               const starling_deriv_options& dopt)
{
    stage_data sdata;

    //
    // HEAD contains everything before the head position in the
    // stage processing pipeline, this should be:
    //
    // 1) GROUPER contigs/reads that were far out-of-order
    // 2) Exon entries that extend past the head position
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

    if (! opt.is_htype_calling)
    {
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

    }
    else
    {
        // POST_REGION defines the end of the "active region" used for
        // haplotyping and calling. the "active region" is contained
        // between two "buffer segments" where haplotyping occurs but
        // no calling is done. The buffer segment is used to improve
        // haplotype continuity on active region boundaries.
        //
        sdata.add_stage(POST_REGION,POST_ALIGN,(opt.htype_buffer_segment()+
                                                opt.htype_call_segment));

        // POST_CALL is used to free up the pileup data. This data is
        // preserved for one extra cycle after snp and indel calling so that
        // the indel caller can look one base 'to the left' of its call
        // location (excepting right breakpoints) to get the depth for the
        // current indel call
        //
        // At the end of this stage we free the position's pileup buffer
        //
        static const unsigned pileup_free_delay(1);
        sdata.add_stage(POST_CALL,POST_REGION,pileup_free_delay);


        // POST_READ is used to free reads during haplotype based calling:
        //
        const unsigned read_reserve_segment(get_read_buffer_size(largest_read_size,
                                                                 largest_total_indel_ref_span_per_read));
        const unsigned post_read_to_post_region(opt.htype_buffer_segment()+read_reserve_segment);
        sdata.add_stage(POST_READ,POST_REGION,post_read_to_post_region);
    }

    // dynamic stages after POST_CALL are used to region based information around a site
    //
    // TODO this will not work correctly for the haplotype calling right now:
    //
    const int pileup_stage(get_pileup_stage_no(opt));

    const auto& stages(dopt.getPostCalLStage());
    const unsigned stageCount(stages.size());
    for (unsigned i(0); i<stageCount; ++i)
    {
        sdata.add_stage(SIZE+i,pileup_stage,stages[i]);
    }

    return sdata;
}
}

starling_pos_processor_base::
starling_pos_processor_base(const starling_options& client_opt,
                            const starling_deriv_options& client_dopt,
                            const reference_contig_segment& ref,
                            const starling_streams_base& client_io,
                            const unsigned n_samples)
    : base_t()
    , _client_opt(client_opt)
    , _client_dopt(client_dopt)
    , _ref(ref)
    , _client_io(client_io)
    , _rmi(STARLING_INIT_LARGEST_READ_SIZE)
    //, _largest_indel_size(std::min(client_opt.max_indel_size,STARLING_INIT_LARGEST_INDEL_SIZE)) -- tmp change for GRUOPER handling
    , _largest_indel_ref_span(client_opt.max_indel_size)
    , _largest_total_indel_ref_span_per_read(_largest_indel_ref_span)
    , _stageman(STAGE::get_stage_data(STARLING_INIT_LARGEST_READ_SIZE, get_largest_total_indel_ref_span_per_read(), _client_opt, _client_dopt),client_dopt.report_range,*this)
    , _chrom_name(_client_opt.bam_seq_name)
    , _n_samples(n_samples)
    , _is_variant_windows(_client_opt.variant_windows.size())
{
    assert((_n_samples != 0) && (_n_samples <= MAX_SAMPLE));

    const unsigned report_size(_client_dopt.report_range.size());
    const unsigned knownref_report_size(get_ref_seq_known_size(_ref,_client_dopt.report_range));
    for (unsigned i(0); i<_n_samples; ++i)
    {
        _sample[i].reset(new sample_info(_client_opt,report_size,knownref_report_size,&_ric));
    }

    // setup gvcf aggregator
    if (_client_opt.is_gvcf_output())
    {
        gvcf_aggregator* temp_agg = new gvcf_aggregator(client_opt,client_dopt,ref,client_io.gvcf_osptr(0));
        if (_client_opt.do_codon_phasing)
        {
            temp_agg->codon_phaser.read_buffer = &_sample[0]->read_buff; // give codon-phaser access to the read-buffer
            temp_agg->codon_phaser.set_options(client_opt,client_dopt);
            temp_agg->codon_phaser.max_read_len = get_largest_read_size();
        }
        _gvcfer.reset(temp_agg);
    }

#ifdef HAVE_FISHER_EXACT_TEST
    if (_client_opt.is_adis_table)
    {
        _ws.reset(get_exact_test_ws());
    }
#endif

    if (_client_opt.is_bsnp_nploid)
    {
        _ninfo.reset(new nploid_info(_client_opt.bsnp_nploid_ploidy));
    }

    if (_client_opt.is_all_sites())
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
            _client_dopt.pdcaller().position_snp_call_pprob_digt(_client_opt,good_epi,
                                                                 *_empty_dgt[b],_client_opt.is_all_sites());
        }
    }

    _is_dependent_eprob = ((_client_opt.is_bsnp_diploid() || _client_opt.is_bsnp_monoploid) &&
                           (_client_opt.bsnp_ssd_no_mismatch>0. || _client_opt.bsnp_ssd_one_mismatch>0));

    // define an expanded indel influence zone around the report range:
    //
    // note that we don't know the max indel ref span per read at this point, so a fudge factor is
    // added here:
    const int bshift(get_influence_zone_size(get_largest_read_size(),
                                             _client_opt.max_indel_size*2));
    pos_range& rir( _report_influence_range);
    rir = _client_dopt.report_range_limit;
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
                              _client_opt,
                              _client_dopt));
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
    assert(is<=_client_opt.max_indel_size);
    _largest_indel_ref_span=std::min(is,_client_opt.max_indel_size);
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

    pos_range output_report_range(_client_dopt.report_range);

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

    if (_client_opt.is_gvcf_output())
    {
        _gvcfer->flush();
    }
}



bool
starling_pos_processor_base::
insert_indel(const indel_observation& obs,
             const unsigned sample_no)
{
    //
    // ppr advance is controlled by the start positions of reads and
    // relatively cheap to store (so long as we aren't includingu
    // gigantic insert sequences) and do not scale up linearly with
    // increased coverage like reads do. For this reason our strategy
    // is to buffer the indels as far ahead as possible while leaving
    // the read buffer sizes fixed at a smaller value.
    //

    try
    {
        _stageman.validate_new_pos_value(obs.key.pos,STAGE::READ_BUFFER);

        // dynamically scale maximum indel size:
        const unsigned len(std::min(static_cast<unsigned>((obs.key.delete_length())),_client_opt.max_indel_size));
        update_largest_indel_ref_span(len);

        bool is_novel(sample(sample_no).indel_sync().insert_indel(obs));
        if (obs.data.is_forced_output) _is_skip_process_pos=false;

        return is_novel;
    }
    catch (...)
    {
        log_os << "Exception caught while attempting to insert indel: " << obs << "\n";
        throw;
    }
}


// snv gt and stats must be reported for this pos (not only honored in strelka right now)
void
starling_pos_processor_base::
insert_forced_output_pos(const pos_t pos)
{
    _stageman.validate_new_pos_value(pos,STAGE::READ_BUFFER);
    _forced_output_pos.insert(pos);
}



bool
starling_pos_processor_base::
is_estimated_depth_range_ge_than(
    const pos_t begin,
    const pos_t end,
    const unsigned depth,
    const unsigned sample_no) const
{
    return sample(sample_no).estdepth_buff.is_range_ge_than(begin,end,depth);
}


// TODO use boost::optional here:
//
std::pair<bool,align_id_t>
starling_pos_processor_base::
insert_read(
    const bam_record& br,
    const alignment& al,
    const char* chrom_name,
    const MAPLEVEL::index_t maplev,
    const unsigned sample_no)
{
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
        return std::make_pair(false,0);
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

    starling_read_buffer& rbuff(sample(sample_no).read_buff);
    const std::pair<bool,align_id_t> res(rbuff.add_read_alignment(_client_opt,
                                                                  br,al,maplev));
    if (! res.first) return res;

    // must initialize initial read_segments "by-hand":
    //
    // TODO get this streamlined into the pos-processor
    //
    {
        const starling_read* sread_ptr(rbuff.get_read(res.second));
        assert(NULL!=sread_ptr);

        // update depth-buffer for the whole read:
        load_read_in_depth_buffer(sread_ptr->get_full_segment(),sample_no);

        // update other data for only the first read segment
        const seg_id_t seg_id(sread_ptr->is_segmented() ? 1 : 0 );
        init_read_segment(sread_ptr->get_segment(seg_id),sample_no);
    }

    return res;
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
        add_alignment_to_depth_buffer(al,sample(sample_no).estdepth_buff);
    }
}

// only acts on genomic mapped reads:
void
starling_pos_processor_base::
init_read_segment(const read_segment& rseg,
                  const unsigned sample_no)
{
    const alignment& al(rseg.genome_align());
    if (al.empty()) return;

    const MAPLEVEL::index_t maplev(rseg.genome_align_maplev());
    //    const bool is_usable_mapping(MAPLEVEL::TIER1_MAPPED == maplev);

    const INDEL_ALIGN_TYPE::index_t iat(translate_maplev_to_indel_type(maplev));

    const bam_seq bseq(rseg.get_bam_read());
    try
    {
//        const read_stats rs = read_stats(rseg.map_qual(),rseg.qual());
        const unsigned total_indel_ref_span_per_read =
            add_alignment_indels_to_sppr(_client_opt.max_indel_size,_ref,
                                         al,bseq,*this,iat,rseg.id(),sample_no,rseg.get_segment_edge_pin());

        update_largest_total_indel_ref_span_per_read(total_indel_ref_span_per_read);

    }
    catch (...)
    {
        log_os << "\nException caught in add_alignment_indels_to_sppr() while processing record: " << rseg.key() << "\n";
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

    const pos_t min_pos(std::max(static_cast<pos_t>(0),pos-static_cast<pos_t>(post_offset-buffer_offset)));
    const pos_t max_pos(pos+1+(buffer_offset-head_offset));
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
    const known_pos_range realign_pr(get_realignment_range(pos, _stageman.get_stage_data()));

    for (unsigned s(0); s<_n_samples; ++s)
    {
        sample_info& sif(sample(s));
        read_segment_iter ri(sif.read_buff.get_pos_read_segment_iter(pos));
        for (read_segment_iter::ret_val r; true; ri.next())
        {
            r=ri.get_ptr();
            if (nullptr == r.first) break;
            read_segment& rseg(r.first->get_segment(r.second));
            if (_client_opt.is_realign_submapped_reads ||
                rseg.is_treated_as_anytier_mapping())
            {
                try
                {
                    realign_and_score_read(_client_opt,_client_dopt,sif.sample_opt,_ref,realign_pr,rseg,sif.indel_sync());
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

        if (_client_opt.is_write_candidate_indels())
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

        if (! _client_opt.is_htype_calling)
        {
            if (! _client_opt.is_write_candidate_indels_only)
            {
                //        clean_pos(pos);
                if (! _client_opt.is_skip_realignment)
                {
                    align_pos(pos);
                }
                pileup_pos_reads(pos);
                // if(_client_opt.is_realigned_read_file) {
                //     rebuffer_pos_reads(pos);
                // }

                write_reads(pos);
            }

            //previously cleared read-buffer here

        }
        else
        {
            if (! _client_opt.is_write_candidate_indels_only)
            {
                //        clean_pos(pos);
                align_pos(pos);
                rebuffer_pos_reads(pos);
            }
        }

    }
    else if (stage_no==STAGE::POST_ALIGN)
    {
        if (! _client_opt.is_htype_calling)
        {
            if (! _client_opt.is_write_candidate_indels_only)
            {
            	if (is_pos_reportable(pos))
                {

                    process_pos_variants(pos);
                }
            }

            // clear read buffer here as oppose to READ_BUFFER stage.
            // if we are doing short-range phasing, the codon_phaser
            // is responsible or clearing the read buffer
            if (!this->_client_opt.do_codon_phasing)
                for (unsigned s(0); s<_n_samples; ++s)
                {
                    sample(s).read_buff.clear_pos(_client_opt.is_ignore_read_names,pos);
                }

            clear_forced_output_pos(pos);

            for (unsigned s(0); s<_n_samples; ++s)
            {
                sample(s).indel_buff.clear_pos(pos);
            }
        }

    }
    else if (stage_no==STAGE::POST_REGION)
    {
        assert(_client_opt.is_htype_calling);

        if (! _client_opt.is_write_candidate_indels_only)
        {
            process_htype_pos(pos);
            clear_forced_output_pos(pos);
        }

    }
    else if (stage_no==STAGE::POST_CALL)
    {
        for (unsigned s(0); s<_n_samples; ++s)
        {
            sample_info& sif(sample(s));
            sif.estdepth_buff.clear_pos(pos);
            sif.bc_buff.clear_pos(pos);
        }
    }
    else if (stage_no==STAGE::POST_READ)
    {
        assert(_client_opt.is_htype_calling);

        for (unsigned s(0); s<_n_samples; ++s)
        {
            sample(s).read_buff.clear_pos(_client_opt.is_ignore_read_names,pos);
        }
        for (unsigned s(0); s<_n_samples; ++s)
        {
            sample(s).indel_buff.clear_pos(pos);
        }

    }
    else if (stage_no>STAGE::POST_CALL)
    {
        run_post_call_step(stage_no,pos);
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
    const uint8_t adjustedMapq,
    const unsigned cycle,
    const bool is_submapped)
{
    // assume pos is already valid:

    auto& bcbuff(sample(sample_no).bc_buff);
    bcbuff.insert_mapq_count(pos,mapq,adjustedMapq);
    bcbuff.update_ranksums(_ref.get_base(pos),pos,call_id,qscore,adjustedMapq,cycle,is_submapped);
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
    bcbuff.insert_mapq_count(pos,mapq,mapq);

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
write_candidate_indels_pos(const pos_t pos)
{
    static const unsigned sample_no(0);

    const pos_t output_pos(pos+1);
    typedef indel_buffer::const_iterator ciiter;

    sample_info& sif(sample(sample_no));
    ciiter i(sif.indel_buff.pos_iter(pos));
    const ciiter i_end(sif.indel_buff.pos_iter(pos+1));

    std::ostream& bos(*_client_io.candidate_indel_osptr());

    for (; i!=i_end; ++i)
    {
        const indel_key& ik(i->first);
        const indel_data& id(get_indel_data(i));
        if (! sif.indel_sync().is_candidate_indel(_client_opt,ik,id)) continue;
        bos << _chrom_name << "\t"
            << output_pos << "\t"
            << INDEL::get_index_label(ik.type) << "\t"
            << ik.length << "\t";
        if (INDEL::SWAP==ik.type)
        {
            bos << ik.swap_dlength << "\t";
        }
        bos << id.get_insert_seq() << "\n";
    }
}



void
starling_pos_processor_base::
process_pos_indel_single_sample(const pos_t pos,
                                const unsigned sample_no)
{
    // note multi-sample status -- can still be called only for one sample
    // and only for sample 0. working on generalization:
    //
    if (sample_no!=0) return;

    // Current multiploid indel model can handle a het or hom indel
    // allele vs. reference, or two intersecting non-reference indel
    // alleles. (note that indel intersection is evaluated only in
    // terms of breakpoints -- so, for instance, a small het deletion
    // could occur within a large het deletion and the two would be
    // treated as non-interacting -- this is just an artifact of how
    // the methods are coded,)
    //

    std::ostream& report_os(std::cerr);
    const pos_t output_pos(pos+1);

    typedef indel_buffer::const_iterator ciiter;

    sample_info& sif(sample(sample_no));
    ciiter it(sif.indel_buff.pos_iter(pos));
    const ciiter it_end(sif.indel_buff.pos_iter(pos+1));

    for (; it!=it_end; ++it)
    {
        const indel_key& ik(it->first);
        const indel_data& id(get_indel_data(it));
        const bool forcedOutput(id.is_forced_output);
        const bool zeroCoverage(id.read_path_lnp.empty());

        if (!sif.indel_sync().is_candidate_indel(_client_opt,ik,id) && !forcedOutput) continue;
        if (zeroCoverage && !forcedOutput) continue;

        // TODO implement indel overlap resolution
        //
        // punt conflict resolution for now....

        if (_client_opt.is_bindel_diploid())
        {
            // indel_report_info needs to be run first now so that
            // local small repeat info is available to the indel
            // caller

            // sample-independent info:
            starling_indel_report_info iri;
            get_starling_indel_report_info(ik,id,_ref,iri);

            double indel_error_prob(0);
            double ref_error_prob(0);
            get_indel_error_prob(_client_opt,iri,indel_error_prob,ref_error_prob);

            static const bool is_tier2_pass(false);
            static const bool is_use_alt_indel(true);

            starling_diploid_indel dindel;
            dindel.is_forced_output = forcedOutput;
            dindel.is_zero_coverage = zeroCoverage;
            _client_dopt.incaller().starling_indel_call_pprob_digt(_client_opt,_client_dopt,
                                                                   sif.sample_opt,
                                                                   indel_error_prob,ref_error_prob,
                                                                   ik,id,is_use_alt_indel,dindel);

            bool is_indel(false);
            if ((dindel.is_indel) or (dindel.is_forced_output))
            {
                is_indel=true;

                // sample-specific info: (division doesn't really matter
                // in single-sample case)
                starling_indel_sample_report_info isri;
                get_starling_indel_sample_report_info(_client_dopt,ik,id,sif.bc_buff,
                                                      is_tier2_pass,is_use_alt_indel,isri);

                if (_client_opt.is_gvcf_output())
                {
                    _gvcfer->add_indel(pos,ik,dindel,iri,isri);
                }

                if (_client_opt.is_bindel_diploid_file)
                {

                    std::ostream& bos(*_client_io.bindel_diploid_osptr(sample_no));
                    bos << _chrom_name << "\t" << output_pos << "\t";
                    write_starling_diploid_indel_file(dindel,iri,isri,bos);
                    bos << "\n";
                }
                if (_is_variant_windows) _variant_print_pos.insert(pos);
            }

            /// \TODO put this option under runtime control...
            /// \TODO setup option so that read keys persist longer when needed for this case...
            ///
            static const bool is_print_indel_evidence(false);

            if (is_print_indel_evidence && is_indel)
            {
                report_os << "INDEL_EVIDENCE " << ik;

                typedef indel_data::score_t::const_iterator siter;
                siter it2(id.read_path_lnp.begin()), it2_end(id.read_path_lnp.end());
                for (; it2!=it2_end; ++it2)
                {
                    const align_id_t read_id(it2->first);
                    const read_path_scores& lnp(it2->second);
                    const read_path_scores pprob(indel_lnp_to_pprob(_client_dopt,lnp,is_tier2_pass,is_use_alt_indel));
                    const starling_read* srptr(sif.read_buff.get_read(read_id));

                    report_os << "read key: ";
                    if (NULL==srptr) report_os << "UNKNOWN_KEY";
                    else            report_os << srptr->key();
                    report_os << "\n"
                              << "read log_lhoods: " << lnp << "\n"
                              << "read pprobs: " << pprob << "\n";
                }
            }
        }
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
        bam_dumper* bamd_ptr(_client_io.realign_bam_ptr(s));
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
    static const bool is_include_submapped(false);

    for (unsigned s(0); s<_n_samples; ++s)
    {
        read_segment_iter ri(sample(s).read_buff.get_pos_read_segment_iter(pos));
        read_segment_iter::ret_val r;
        while (true)
        {
            r=ri.get_ptr();
            if (nullptr==r.first) break;
            const read_segment& rseg(r.first->get_segment(r.second));
            if (is_include_submapped || rseg.is_treated_as_anytier_mapping())
            {
                pileup_read_segment(rseg,s);
            }
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
        if (! rseg.is_any_nonovermax(_client_opt.max_indel_size)) return;
    }

    const alignment& best_al(*best_al_ptr);
    if (best_al.empty())
    {
        if (! rseg.is_realigned)
        {
            if (_client_opt.verbosity >= LOG_LEVEL::ALLWARN)
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
    const bool is_mapq_adjust(adjustedMapq<=80);
    // test read against max indel size (this is a backup, should have been taken care of upstream):
    const unsigned read_ref_mapped_size(apath_ref_length(best_al.path));
    if (read_ref_mapped_size > (read_size+get_largest_total_indel_ref_span_per_read()))
    {
        //brc.large_ref_deletion++;
        return;
    }

    // exact begin and end report range filters:
    {
        const pos_range& rlimit(_client_dopt.report_range_limit);
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

    const bool is_submapped(! rseg.is_treated_as_anytier_mapping());
    const bool is_tier1(rseg.is_treated_as_tier1_mapping());

    // precompute mismatch density info for this read:
    if ((! is_submapped) && _client_opt.is_max_win_mismatch)
    {
        const rc_segment_bam_seq ref_bseq(_ref);
        create_mismatch_filter_map(_client_opt,best_al,ref_bseq,bseq,read_begin,read_end,_rmi);
        if (_client_opt.is_tier2_mismatch_density_filter_count)
        {
            const int max_pass(_client_opt.tier2_mismatch_density_filter_count);
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

    using namespace ALIGNPATH;
    const unsigned as(best_al.path.size());
    for (unsigned i(0); i<as; ++i)
    {
        const path_segment& ps(best_al.path[i]);

#ifdef DEBUG_PPOS
        log_os << "seg,ref,read: " << i << " " << ref_head_pos << " " << read_head_pos << "\n";
#endif

        if (ps.type == MATCH)
        {
            for (unsigned j(0); j<ps.length; ++j)
            {
                const unsigned read_pos(read_head_pos+j);
                if ((read_pos < read_begin) || (read_pos >= read_end)) continue; // allow for read end trimming
                const pos_t ref_pos(ref_head_pos+static_cast<pos_t>(j));

#ifdef DEBUG_PPOS
                log_os << "j,ref,read: " << j << " " << ref_pos <<  " " << read_pos << "\n";
#endif

                // skip position outside of report range:
                if (! is_pos_reportable(ref_pos)) continue;
                _stageman.validate_new_pos_value(ref_pos,STAGE::get_pileup_stage_no(_client_opt));

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
                                        (qscore < _client_opt.min_qscore));

                    assert(! _client_opt.is_min_win_qscore);

                    bool is_tier2_call_filter(is_call_filter);
                    if (! is_call_filter && _client_opt.is_max_win_mismatch)
                    {
                        is_call_filter = _rmi[read_pos].mismatch_filter_map;
                        if (! _client_opt.is_tier2_no_mismatch_density_filter)
                        {
                            if (_client_opt.is_tier2_mismatch_density_filter_count)
                            {
                                is_tier2_call_filter = _rmi[read_pos].tier2_mismatch_filter_map;
                            }
                            else
                            {
                                is_tier2_call_filter = is_call_filter;
                            }
                        }
                    }
                    current_call_filter = ( is_tier1 ? is_call_filter : is_tier2_call_filter );
                    is_tier_specific_filter = ( is_tier1 && is_call_filter && (! is_tier2_call_filter) );

                    if (_client_opt.is_max_win_mismatch)
                    {
                        is_neighbor_mismatch=(_rmi[read_pos].mismatch_count_ns>0);
                    }
                }

                // update extended feature metrics (including submapped reads):
                if (_client_opt.is_compute_germline_VQSRmetrics())
                {
                    /// \TODO Morten -- consider improving MQ, MQ0 and RankSumMQ by:
                    ///  1) removing the if (! submapped) here
                    ///  2) collapsing mapq and adjustedMapq to just mapq
                    ///
                    if (! is_submapped)
                    {
                        update_ranksum_and_mapq_count(ref_pos,sample_no,call_id,qscore,mapq,adjustedMapq,align_strand_read_pos,is_submapped);
                    }
                }
                else if (_client_opt.is_compute_somatic_VQSRmetrics)
                {
                    update_somatic_features(ref_pos,sample_no,is_tier1,call_id,current_call_filter,mapq,read_pos,read_size);
                }

                if (is_submapped)
                {
                    insert_pos_submap_count(ref_pos,sample_no);
                    continue;
                }

                /// include only data meeting mapping criteria after this point:

                if (_client_opt.is_compute_hapscore)
                {
                    insert_hap_cand(ref_pos,sample_no,is_tier1,
                                    bseq,qual,read_pos);
                }

                try
                {
                    const base_call bc = base_call(call_id,qscore,best_al.is_fwd_strand,
                                                   align_strand_read_pos,end_trimmed_read_len,
                                                   current_call_filter,is_neighbor_mismatch,is_tier_specific_filter);

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
                    _stageman.validate_new_pos_value(ref_pos,STAGE::get_pileup_stage_no(_client_opt));

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
process_pos_snp_single_sample(const pos_t pos,
                              const unsigned sample_no)
{
    try
    {
        process_pos_snp_single_sample_impl(pos,sample_no);
    }
    catch (...)
    {
        log_os << "Exception caught in starling_pos_processor_base.process_pos_snp_single_sample_impl() while processing chromosome position: " << (pos+1) << "\n"
               << "snp_pos_info:\n";
        const snp_pos_info* spi_ptr(sample(sample_no).bc_buff.get_pos(pos));
        if (NULL==spi_ptr)
        {
            static const snp_pos_info spi_null;
            spi_ptr=&spi_null;
        }
        log_os << *spi_ptr << "\n";
        throw;
    }
}



void
starling_pos_processor_base::
process_pos_snp_single_sample_impl(
    const pos_t pos,
    const unsigned sample_no)
{
    // TODO:
    //
    // note this might not matter wrt larger changes taking place, but here goes:
    //
    // change filters to support vcf concept of 1..N filters which are added to the genotype information
    //
    // generalize site tests with an object
    //
    // genotype_test {
    //    ctor(); // setup any cached values
    //
    //    test(site_info);
    //
    //    write()?? (do we need to even bother with this?)
    // }
    //

    sample_info& sif(sample(sample_no));

    snp_pos_info null_pi;
    snp_pos_info* pi_ptr(sif.bc_buff.get_pos(pos));
    if (NULL==pi_ptr) pi_ptr=&null_pi;
    snp_pos_info& pi(*pi_ptr);

    const unsigned n_calls(pi.calls.size());
    const unsigned n_spandel(pi.n_spandel);
    const unsigned n_submapped(pi.n_submapped);

    const pos_t output_pos(pos+1);

    pi.set_ref_base(_ref.get_base(pos));

    // for all but coverage-tests, we use a high-quality subset of the basecalls:
    //
    snp_pos_info& good_pi(sif.epd.good_pi);
    good_pi.clear();
    good_pi.set_ref_base(pi.get_ref_base());
    for (unsigned i(0); i<n_calls; ++i)
    {
        if (pi.calls[i].is_call_filter) continue;
        good_pi.calls.push_back(pi.calls[i]);
    }

    _site_info.n_used_calls=(good_pi.calls.size());
    _site_info.n_unused_calls=(n_calls-_site_info.n_used_calls);

    sif.ss.update(n_calls);
    sif.used_ss.update(_site_info.n_used_calls);
    if (pi.get_ref_base() != 'N')
    {
        sif.ssn.update(n_calls);
        sif.used_ssn.update(_site_info.n_used_calls);
        sif.wav.insert(pos,_site_info.n_used_calls,_site_info.n_unused_calls,n_spandel,n_submapped);
    }
    else
    {
        sif.wav.insert_null(pos);
    }

    // note multi-sample status -- can still be called only for one sample
    // and only for sample 0. working on generalization:
    //
    if (sample_no!=0) return;

    if (pi.calls.empty()) return;

    adjust_joint_eprob(_client_opt,_dpcache,good_pi,sif.epd.dependent_eprob,_is_dependent_eprob);

    const extended_pos_info good_epi(good_pi,sif.epd.dependent_eprob);

    // get fraction of filtered bases:
#if 0
    const double filter_fraction(static_cast<double>(n_unused_calls)/static_cast<double>(n_calls));
    const bool is_overfilter(filter_fraction > _client_opt.max_basecall_filter_fraction);
#endif

    // delay writing any snpcalls so that anomaly tests can (optionally) be applied as filters:
    //
    nonref_test_call nrc;
    //lrt_snp_call lsc;
    _site_info.dgt.reset();
    //monoploid_genotype mgt;
    //std::unique_ptr<nploid_genotype> ngt_ptr;

    if (_client_opt.is_counts)
    {
        report_counts(good_pi,_site_info.n_unused_calls,output_pos,*_client_io.counts_osptr());
    }

    if (_client_opt.is_nonref_test() || _client_opt.is_nonref_sites())
    {
        position_nonref_2allele_test(good_pi,
                                     _client_opt,
                                     _client_opt.is_nonref_sites(),
                                     nrc);
#if 0
        static const bool is_mle_freq(false);

        position_nonref_test(good_pi,
                             _client_opt.nonref_variant_rate,
                             _client_opt.min_nonref_freq,
                             is_mle_freq,
                             nrc);
#endif

    }

#if 0
    if (_client_opt.is_lsnp)
    {
        position_snp_call_lrt(_client_opt.lsnp_alpha,good_pi,lsc);
    }
#endif
    if (_client_opt.is_bsnp_diploid())
    {
        _client_dopt.pdcaller().position_snp_call_pprob_digt(_client_opt,good_epi,_site_info.dgt,_client_opt.is_all_sites());
    }
#if 0
    if (_client_opt.is_bsnp_monoploid)
    {
        position_snp_call_pprob_monogt(_client_opt.bsnp_monoploid_theta,good_pi,mgt);
    }
    if (_client_opt.is_bsnp_nploid)
    {
        ngt_ptr.reset(new nploid_genotype(*_ninfo));
        position_snp_call_pprob_nploid(_client_opt.bsnp_nploid_snp_prob,good_pi,*_ninfo,*ngt_ptr);
    }
#endif

    //    const bool is_snp(nrc.is_snp || lsc.is_snp || _site_info.dgt.is_snp || mgt.is_snp || (ngt_ptr.get() && ngt_ptr->is_snp));
    const bool is_snp(nrc.is_snp || _site_info.dgt.is_snp);

    // find anomalies:
    //
#if 0
    bool is_pos_adis(false);
    bool is_pos_acov(false);

    if ((_client_opt.is_adis_table || _client_opt.is_adis_lrt) && is_snp)
    {
        if (_client_opt.is_adis_table)
        {
            is_pos_adis = (is_pos_adis || position_strand_distro_anomaly(_client_opt.adis_table_alpha,good_pi,_ws));
        }
        if (_client_opt.is_adis_lrt)
        {
            is_pos_adis = (is_pos_adis || position_strand_distro_anomaly_lrt(_client_opt.adis_lrt_alpha,good_pi));
        }
    }
    if (_client_opt.is_acov)
    {
        is_pos_acov = position_strand_coverage_anomaly(_client_opt.acov_alpha,pi);
    }
#endif

    //const bool is_anomaly(is_pos_adis || is_pos_acov);
    //const bool is_filter_snp(is_overfilter || (_client_opt.is_filter_anom_calls && is_anomaly));

    //    const bool is_nf_snp(is_snp && (! is_filter_snp));
    if (is_snp)
    {
        if (_client_opt.is_compute_hapscore)
        {
            _site_info.hapscore=get_hapscore(pi.hap_set);
        }

        // do calculate VQSR metrics
        if (_client_opt.is_compute_germline_VQSRmetrics())
        {
            _site_info.MQ 				= pi.get_rms_mq();
            _site_info.ReadPosRankSum 	= pi.get_read_pos_ranksum();
            _site_info.MQRankSum 		= pi.get_mq_ranksum();
            _site_info.BaseQRankSum 	= pi.get_baseq_ranksum();
            _site_info.rawPos           = pi.get_raw_pos();
            _site_info.avgBaseQ         = pi.get_raw_baseQ();
        }

        // hpol filter
        _site_info.hpol=get_snp_hpol_size(pos,_ref);
    }

    if (_client_opt.is_all_sites())
    {
#if 0
        const diploid_genotype* dgt_ptr(&_site_info.dgt);
        if (is_filter_snp)
        {
            dgt_ptr=&get_empty_dgt(pi.ref_base);
        }
#endif

        //Add site to gvcf
        if (_client_opt.is_gvcf_output())
        {
            _site_info.init(pos,pi.get_ref_base(),good_pi,_client_opt.used_allele_count_min_qscore);
            _gvcfer->add_site(_site_info);
        }


        if (_client_opt.is_bsnp_diploid_allele_file)
        {
            write_bsnp_diploid_allele(_client_opt,_client_io,_chrom_name,output_pos,pi.get_ref_base(),_site_info.n_used_calls,_site_info.n_unused_calls,good_pi,_site_info.dgt,_site_info.hpol);
        }
    }

    if (_client_opt.is_nonref_sites())
    {
        std::ostream& bos(*_client_io.nonref_sites_osptr());
        write_snp_prefix_info_file(_chrom_name,output_pos,pi.get_ref_base(),_site_info.n_used_calls,_site_info.n_unused_calls,bos);
        bos << "\t";
        write_nonref_2allele_test(_client_opt,good_pi,nrc,bos);
        bos << "\n";
    }

    // report events:
    //
    bool is_reported_event(false);

    std::ostream& report_os(std::cerr);

    if (is_snp)
    {
        if (nrc.is_snp)
        {
            std::ostream& bos(*_client_io.nonref_test_osptr());
            write_snp_prefix_info_file(_chrom_name,output_pos,pi.get_ref_base(),_site_info.n_used_calls,_site_info.n_unused_calls,bos);
            bos << "\t";
            write_nonref_2allele_test(_client_opt,good_pi,nrc,bos);
#if 0
            write_nonref_test(_client_opt,good_pi,nrc,bos);
#endif
            bos << "\n";
        }
#if 0
        if (lsc.is_snp)
        {
            write_snp_prefix_info("LSNP",output_pos,pi.ref_base,_site_info.n_used_calls,_site_info.n_unused_calls,report_os);
            report_os << " " << lsc << "\n";
        }
#endif
        if (_site_info.dgt.is_snp)
        {
            if (_client_opt.is_bsnp_diploid_file)
            {
                std::ostream& bos(*_client_io.bsnp_diploid_osptr());
                write_snp_prefix_info_file(_chrom_name,output_pos,pi.get_ref_base(),_site_info.n_used_calls,_site_info.n_unused_calls,bos);
                bos << "\t";
                write_diploid_genotype_snp(_client_opt,good_pi,_site_info.dgt,bos,_site_info.hpol);
                bos << "\n";
            }

            // this needs to be updated no matter where the snp-call is written to:
            if (_is_variant_windows) _variant_print_pos.insert(pos);
        }
#if 0
        if (mgt.is_snp)
        {
            write_snp_prefix_info("BSNP1",output_pos,pi.ref_base,_site_info.n_used_calls,_site_info.n_unused_calls,report_os);
            report_os << " " << mgt << "\n";
        }
        if (ngt_ptr.get() && ngt_ptr->is_snp)
        {
            write_snp_prefix_info("BSNPN",output_pos,pi.ref_base,_site_info.n_used_calls,_site_info.n_unused_calls,report_os);
            report_os << " ";
            nploid_write(*_ninfo,*ngt_ptr,report_os);
            report_os << "\n";
        }
#endif

        is_reported_event = true;
    }

#if 0
    if (is_anomaly && (! _client_opt.is_filter_anom_calls))
    {
        if (is_pos_adis) report_os << "ANOM_DIS pos: " << output_pos << "\n";
        if (is_pos_acov) report_os << "ANOM_COV pos: " << output_pos << "\n";

        is_reported_event = true;
    }
#endif

    if (_client_opt.is_print_all_site_evidence || (_client_opt.is_print_evidence && is_reported_event))
    {
        report_os << "EVIDENCE pos: " << output_pos << "\n"
                  << "is_snp: " << is_snp << "\n"
                  << pi << "\n";
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



void
starling_pos_processor_base::
run_post_call_step(
    const int stage_no,
    const pos_t pos)
{
    if (_variant_print_pos.count(pos)==0) return;

    const int pcn(STAGE::get_last_static_stage_no(_client_opt));
    assert(stage_no>pcn);

    // convert stage_no to window_no:
    assert(stage_no >= static_cast<int>(_client_dopt.variant_window_first_stage));
    assert(stage_no <= static_cast<int>(_client_dopt.variant_window_last_stage));

    const unsigned window_no(stage_no-_client_dopt.variant_window_first_stage);
    const unsigned vs(_client_opt.variant_windows.size());

    assert(window_no<vs);

    const pos_t output_pos(pos+1);
    std::ostream& bos(*_client_io.variant_window_osptr(window_no));

    bos << _chrom_name << '\t' << output_pos;

    std::ofstream tmp_os;
    tmp_os.copyfmt(bos);

    bos << std::setprecision(2) << std::fixed;

    for (unsigned s(0); s<_n_samples; ++s)
    {
        const win_avg_set& was(sample(s).wav.get_win_avg_set(window_no));
        bos << '\t' << was.ss_used_win.avg()
            << '\t' << was.ss_filt_win.avg()
            << '\t' << was.ss_submap_win.avg();
    }

    bos.copyfmt(tmp_os);

    bos << '\n';

    if ((window_no+1) == vs)
    {
        _variant_print_pos.erase(pos);
    }
}


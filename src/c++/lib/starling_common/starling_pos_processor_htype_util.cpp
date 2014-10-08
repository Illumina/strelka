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

/// \file
///
/// part of starling_pos_processor_base that handles new haplotype routines
///
/// \author Chris Saunders
///

#include "alignment_util.hh"
#include "blt_util/log.hh"
#include "starling_common/starling_indel_report_info.hh"
#include "starling_common/starling_pos_processor_base.hh"

#include <cstdlib>

#include <iostream>



// used to describe an open breakpoint haplotype:
//
namespace OPEN
{
enum index_t
{
    NONE,
    LEFT,
    RIGHT
};

static
const char*
label(index_t i)
{
    switch (i)
    {
    case LEFT:
        return "left";
    case RIGHT:
        return "right";
    default:
        return "";
    }
}
}

namespace
{

// haplotypes are composed of sequences of "elements", each element is
// similar to the current "sequence swap" indel type, thus it can
// represent a snp:
//
// Note that one major difference between the indels and the htype
// element is that any inserted sequence is part of the type and
// otherwise identical elements with different sequences will not
// compare as equal.
//
// There should be no need to have reference elements, as these are
// reflected in a haplotype simply by having all overlapping elements
// turned-off. The presense of an element always implies that we will
// test the corresponding reference haplotype.
//
// All conventions for pos,right_pos(), etc should match indel_key
// when possible
//
struct htype_element
{

    htype_element(const pos_t p=0)
        : pos(p)
        , delete_length(0)
        , open_end(OPEN::NONE) {}

    void
    clear()
    {
        pos=0;
        delete_length=0;
        open_end=OPEN::NONE;
        seq.clear();
    }

    unsigned insert_length() const
    {
        return seq.size();
    }

    pos_t right_pos() const
    {
        return (pos+delete_length);
    }

    bool
    operator==(const htype_element& rhs) const
    {
        return ((pos == rhs.pos) &&
                (delete_length == rhs.delete_length) &&
                (open_end == rhs.open_end) &&
                (seq == rhs.seq));
    }

    bool
    operator<(const htype_element& rhs) const
    {
        if (pos < rhs.pos) return true;
        if (pos==rhs.pos)
        {
            return gtcore(rhs);
        }
        return false;
    }

    bool
    gtcore(const htype_element& rhs) const
    {
        if (delete_length < rhs.delete_length) return true;
        if (delete_length==rhs.delete_length)
        {
            if (open_end < rhs.open_end) return true;
            if (open_end==rhs.open_end)
            {
                if (seq < rhs.seq) return true;
            }
        }
        return false;
    }

    // correct pos range to use when we view sv's as breakpoints:
    known_pos_range
    breakpoint_pos_range() const
    {
        return known_pos_range(pos,right_pos());
    }

    // correct pos range to use when we view sv's as ranges
    // (ie. candidate indel interference within a single read:)
    pos_range
    open_pos_range() const
    {
        if       (open_end == OPEN::RIGHT)
        {
            pos_range pr;
            pr.set_begin_pos(pos);
            return pr;
        }
        else if (open_end == OPEN::LEFT)
        {
            pos_range pr;
            pr.set_end_pos(pos);
            return pr;
        }

        return breakpoint_pos_range();
    }

    bool is_breakpoint() const
    {
        return (open_end != OPEN::NONE);
    }

    pos_t pos;
    unsigned delete_length;
    OPEN::index_t open_end;
    std::string seq; // TODO -- turn this into a bam_seq...
};


}


namespace
{

struct right_pos_htype_element_sorter
{
    bool
    operator()(const htype_element& h1,
               const htype_element& h2) const
    {
        if (h1.right_pos() < h2.right_pos()) return true;
        if (h1.right_pos() == h2.right_pos())
        {
            return h1.gtcore(h2);
        }
        return false;
    }
};

}

std::ostream& operator<<(std::ostream& os, const htype_element& he);


std::ostream&
operator<<(std::ostream& os,
           const htype_element& he)
{

    os << "htype_element: pos: " << he.pos
       << " del_len: " << he.delete_length;
    if (OPEN::NONE != he.open_end)
    {
        os << " open_end: " << OPEN::label(he.open_end);
    }
    os << " seq: " << he.seq;

    return os;
}



namespace
{

struct htype_buffer
{
    typedef std::map<htype_element,unsigned> hdata_t;
    typedef hdata_t::iterator iterator;
    typedef hdata_t::const_iterator const_iterator;

    bool
    empty() const
    {
        return _hdata.empty();
    }

    //
    //
    void
    insert_element(const htype_element& he)
    {
        const hdata_t::iterator i(_hdata.find(he));
        if (i==_hdata.end())
        {
            _hdata.insert(std::make_pair(he,1));
            _rightpos.insert(std::make_pair(he.right_pos(),he.pos));
        }
        else
        {
            i->second++;
        }
    }

    // position iterators which return (at least) all haplotypes with
    // a left or right breakpoint in the range.
    //
    // Note:
    // 1) indels which encompass the range are not returned
    // 2) some non-intersecting indels may be returned in the
    //    iteration range
    //
    std::pair<iterator,iterator>
    pos_range_iter(const pos_t begin_pos, const pos_t end_pos);

    void
    dump(std::ostream& os) const;

private:
#if 0
    typedef boost::multi_index_container<htype_element,
            indexed_by<
            ordered_unique<indentity<htype_element> >,
            ordered_unique<tag<right_pos>,indentity<htype_element>, right_pos_htype_element_sorter>
            >
            > hdata_t;
#endif

    pos_t
    leftmost_rightkey_pos(const pos_t& begin_range_pos,
                          const pos_t& end_range_pos) const;

    hdata_t _hdata;

    typedef std::map<pos_t,pos_t> rightkey_t;
    rightkey_t _rightpos;
};

}


// search for a key in the right_pos-sorted list with value less
// than begin_key:
pos_t
htype_buffer::
leftmost_rightkey_pos(const pos_t& begin_range_pos,
                      const pos_t& end_range_pos) const
{

    typedef rightkey_t::const_iterator riter;

    pos_t begin_pos(end_range_pos);
    riter ri(_rightpos.lower_bound(begin_range_pos));
    const riter ri_end(_rightpos.lower_bound(end_range_pos));
    for (; ri!=ri_end; ++ri) begin_pos=std::min(begin_pos,ri->second);
    return begin_pos;
}



std::pair<htype_buffer::iterator,htype_buffer::iterator>
htype_buffer::
pos_range_iter(const pos_t begin_pos,
               const pos_t end_pos)
{

    typedef htype_element hkey_t;

    const hkey_t end_range_key(end_pos);
    const iterator end(_hdata.lower_bound(end_range_key));
    const hkey_t begin_range_key(begin_pos);
    iterator begin(_hdata.lower_bound(begin_range_key));
    hkey_t left_begin_key(end_range_key);
    if ((begin!=_hdata.end()) && (begin->first<left_begin_key)) left_begin_key=begin->first;
    const hkey_t begin_key(leftmost_rightkey_pos(begin_range_key.pos,end_range_key.pos));
    if (begin_key < left_begin_key) begin=_hdata.find(begin_key);
    return std::make_pair(begin,end);
}



void
htype_buffer::
dump(std::ostream& os) const
{

    os << "Haplotype buffer dump ON\n";
    const_iterator i(_hdata.begin()),i_end(_hdata.end());
    for (; i!=i_end; ++i)
    {
        os << i->first << " count: " << i->second << "\n";
    }
    os << "Haplotype buffer dump OFF\n";
}



// 99% of this task is taking care of indel normalization
static
bool
convert_indel_to_htype(const indel_key& ik,
                       const indel_data& /*id*/,
                       const read_segment& rseg,
                       const reference_contig_segment& ref,
                       htype_element& he)
{

    he.clear();

    // get best alignment:
    const alignment* alptr(rseg.get_best_alignment());

    assert(alptr);
    const alignment& al(*alptr);

    // Check that alignment is compatible with indel. Many
    // cases where this fails will be for 'private'
    // indels. The posterior above is over all candidate
    // indels, so one candidate may be the best for this read,
    // *but* the best alignment contains a private indel
    // instead.
    //
    pos_range read_indel_pr;
    if (! is_indel_in_alignment(al,ik,read_indel_pr)) return false;

    const bam_seq read_seq(rseg.get_bam_read());

    const rc_segment_bam_seq ref_bseq(ref);
    pos_range ref_indel_pr(ik.open_pos_range());

    assert(! read_indel_pr.is_empty());
    assert(! ref_indel_pr.is_empty());

    // normalization function adjusts ranges:
    //    normalize_indel(read_indel_pr,ref_indel_pr,read_seq,ref_bseq,read_indel_pr,ref_indel_pr);

    assert(! read_indel_pr.is_empty());
    assert(! ref_indel_pr.is_empty());

    // build he:
    if (ref_indel_pr.is_complete())
    {
        he.delete_length=ref_indel_pr.end_pos-ref_indel_pr.begin_pos;
    }

    if (!  read_indel_pr.is_begin_pos)
    {
        he.pos=read_indel_pr.end_pos;
    }
    else
    {
        he.pos=read_indel_pr.begin_pos;
    }

    if       (! read_indel_pr.is_end_pos)
    {
        he.open_end=OPEN::RIGHT;
    }
    else if (! read_indel_pr.is_begin_pos)
    {
        he.open_end=OPEN::LEFT;
    }

    {
        // copy into htype element seq (don't worry about efficiency for now)
        pos_range pr(read_indel_pr);
        if (! pr.is_complete())
        {
            const pos_range nonclip_pr(get_nonclip_range(al.path));
            assert(nonclip_pr.is_complete());
            if       (! pr.is_begin_pos)
            {
                pr.set_begin_pos(nonclip_pr.begin_pos);
            }
            else
            {
                pr.set_end_pos(nonclip_pr.end_pos);
            }
        }

        assert(pr.begin_pos<=pr.end_pos && pr.begin_pos>=0);
        for (pos_t i(pr.begin_pos); i<pr.end_pos; ++i)
        {
            he.seq.push_back(read_seq.get_char(i));
        }
    }

    if ((he.delete_length==0) &&
        (he.insert_length()==0))
    {
        he.clear();
        return false;
    }

    return true;
}




static
void
get_htypes_for_indel(const starling_deriv_options& dopt,
                     const starling_read_buffer& rbuff,
                     const reference_contig_segment& ref,
                     const indel_key& ik,
                     const indel_data& id,
                     htype_buffer& hdata)
{

    static const bool is_tier2_pass(false);
    static const bool is_use_alt_indel(true);
    static const double include_thresh(0.999);

    typedef indel_data::score_t::const_iterator siter;

    htype_element he;
    siter i(id.read_path_lnp.begin()), i_end(id.read_path_lnp.end());
    for (; i!=i_end; ++i)
    {
        const align_id_t read_id(i->first);
        const read_path_scores& path_lnp(i->second);
        const read_path_scores pprob(indel_lnp_to_pprob(dopt,path_lnp,is_tier2_pass,is_use_alt_indel));

        if (pprob.indel >= include_thresh)
        {
            // convert indel to htype_element and insert:
            const starling_read* srptr(rbuff.get_read(read_id));

            /// TODO - add segmented read support
            if (srptr->is_segmented())
            {
                log_os << "ERROR: haplotype model does not work for spliced reads\n";
                exit(EXIT_FAILURE);
            }
            const read_segment& rseg(srptr->get_full_segment());

            if (! convert_indel_to_htype(ik,id,rseg,ref,he)) continue;

            hdata.insert_element(he);
        }
    }
}



// initial step:
// 0) id all indels in full_pr
// 1) convert these into genotype elements (including simple homopolymer normalization)
// 2) reduce to two most popular alleles for all alleles which conflict
// 3) print lots of stuff out for the active region
//
void
starling_pos_processor_base::
get_region_haplotypes(const known_pos_range full_pr,
                      const known_pos_range /*active_pr*/)
{

    // only works for the first sample right now:
    static const unsigned sample_no(0);

    //for(unsigned sample_no(0);sample_no<_n_samples;++sample_no)..
    sample_info& sif(sample(sample_no));

    // step 0: identify all candidate indels in full region:
    htype_buffer hdata;

    // iterate through indels in each region, find and iterate through supporting reads for each indel:
    //
    const indel_buffer& ibuff(sif.indel_sync().ibuff());
    typedef indel_buffer::const_iterator ciiter;
    const std::pair<ciiter,ciiter> ipair(ibuff.pos_range_iter(full_pr.begin_pos,full_pr.end_pos));
    for (ciiter i(ipair.first); i!=ipair.second; ++i)
    {
        const indel_key& ik(i->first);
        const indel_data& id(get_indel_data(i));

        // check that indel really falls into candidate range:
        if (! (full_pr.is_pos_intersect(ik.pos) ||
               full_pr.is_pos_intersect(ik.right_pos()))) continue;

        // only consider candidates:
        if (! sif.indel_sync().is_candidate_indel(ik,id)) continue;

        get_htypes_for_indel(_client_dopt,sif.read_buff,_ref,ik,id,hdata);

        // convert indel to htype_element and insert:
        //  convert_indel_to_htype(ik,id,_ref_seq);
    }

#ifdef DEBUG_HTYPE
    if (! hdata.empty()) hdata.dump(std::cerr);
#endif

#if 0
    // get expanded region for read search:
    const known_pos_range expanded_pr(full_pr.begin_pos-(_client_opt.max_indel_size+XXXMAX_READ_SIZE),full_pr.end_pos);

    // iterate through reads in region:
    for (pos_t pos(begin); pos<end; ++pos)
    {
        read_segment_iter ri(sif.read_buff.get_pos_read_segment_iter(pos));
        for (read_segment_iter::ret_val r; true; ri.next())
        {
            r=ri.get_ptr();
            if (NULL==r.first) break;
            read_segment& rseg(r.first->get_segment(r.second));

            {
                indel_set_t cal_indels;
                get_alignment_indels(cal,opt.max_indel_size,cal_indels);
            }

            //rseg.blah....
        }
    }
#endif

#if 0
    const indel_buffer& ibuff(sif.indel_sync().ibuff());
    const std::pair<ciiter,ciiter> ipair(ibuff.pos_range_iter(full_pr.begin_pos,full_pr.end_pos));
    for (ciiter i(ipair.first); i!=ipair.second; ++i)
    {
        const indel_key& ik(i->first);
        const indel_data& id(get_indel_data(i));

        // check that indel really falls into candidate range:
        if (! (full_pr.is_pos_intersect(ik.pos()) ||
               full_pr.is_pos_intersect(ik.right_pos()))) continue;

        // only consider candidates:
        if (! sif.indel_sync().is_candidate_indel(_client_opt,ik,id)) continue;

        // convert indel to htype_element and insert:
        convert_indel_to_htype(ik,id,_ref_seq);
    }
#endif


#if 0
    const known_pos_range full_pr(start_pos,start_pos+region_size);
    const known_pos_range active_pr(full_start_pos+region_size/4,start_pos+region_size);
    const pos_t begin(region_size/4);
    const pos_t end(begin+(region_size/2));

    // step 1: identify loci in region:
    const indel_buffer& ibuff(isync.ibuff());
    typedef indel_buffer::const_iterator ciiter;
    const std::pair<ciiter,ciiter> ipair(ibuff.pos_range_iter(pr.begin_pos,pr.end_pos));
    for (ciiter i(ipair.first); i!=ipair.second; ++i)
    {
        const indel_key& ik(i->first);
        const indel_data& id(get_indel_data(i));
#if 1
        // iterate through reads in region:
        for (pos_t pos(begin); pos<end; ++pos)
        {
            read_segment_iter ri(sif.read_buff.get_pos_read_segment_iter(pos));
            for (read_segment_iter::ret_val r; true; ri.next())
            {
                r=ri.get_ptr();
                if (NULL==r.first) break;
                read_segment& rseg(r.first->get_segment(r.second));

                //rseg.blah....
            }
        }
#endif
    }
#endif
}



void
starling_pos_processor_base::
process_htype_pos(const pos_t begin_pos)
{

    // exclude negative begin_pos b/c
    // (1) don't wan't to have to write negative mod
    // (2) should be disallowed in any BAM-based pipeline
    if (begin_pos<0) return;

    if (_hregion.is_first_region)
    {
        _hregion.is_first_region=false;
        _hregion.region_alignment=begin_pos%static_cast<pos_t>(_client_opt.htype_call_segment);
    }
    else
    {
        if ((begin_pos%static_cast<pos_t>(_client_opt.htype_call_segment))!=_hregion.region_alignment) return;
    }

    const known_pos_range active_pr(begin_pos,begin_pos+_client_opt.htype_call_segment);

    // see if there's anything reportable in the central region:
    bool is_reportable(false);

    for (pos_t pos(active_pr.begin_pos); pos<active_pr.end_pos; ++pos)
    {
        if (is_pos_reportable(pos))
        {
            is_reportable=true;
            break;
        }
    }
    if (! is_reportable) return;

    // do haplotypeing routines:
    const unsigned full_htype_segment(_client_opt.htype_buffer_segment()+
                                      _client_opt.htype_call_segment+
                                      _client_opt.htype_buffer_segment());

    const pos_t full_begin_pos(begin_pos-_client_opt.htype_buffer_segment());
    const known_pos_range full_pr(full_begin_pos,full_begin_pos+full_htype_segment);
    get_region_haplotypes(full_pr,active_pr);

    // do regular calling for now:
    for (pos_t pos(active_pr.begin_pos); pos<active_pr.end_pos; ++pos)
    {
        pileup_pos_reads(pos);
        write_reads(pos);
        if (is_pos_reportable(pos))
        {
            process_pos_variants(pos);
        }
    }
}


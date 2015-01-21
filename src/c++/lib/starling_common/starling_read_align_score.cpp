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
/// \author Chris Saunders
///

#include "starling_read_align_score.hh"

#include "blt_util/blt_exception.hh"
#include "blt_util/log.hh"
#include "blt_util/qscore.hh"
#include "blt_util/seq_util.hh"
#include "blt_util/align_path_util.hh"

#include <cassert>

#include <sstream>



//#define DEBUG_SCORE



#ifdef DEBUG_SCORE
#include <sstream>
#include <iostream>

struct align_position
{

    align_position(const char r,
                   const char f,
                   const char i) : read(r), ref(f), insert(i) {}
    char read;
    char ref;
    char insert;
};


struct align_printer
{

    void
    push(const char r, const char f, const char i)
    {
        _seq.push_back(align_position(r,f,i));
    }

    void
    dump(std::ostream& os) const
    {
        os << "scoring alignment:\n";
        const unsigned ss(_seq.size());
        os << "read:   ";
        for (unsigned i(0); i<ss; ++i) os << _seq[i].read;
        os << "\n";
        os << "        ";
        for (unsigned i(0); i<ss; ++i)
        {
            char c(' ');
            if (_seq[i].read != '-')
            {
                char r(_seq[i].ref);
                if (r == '-') r=_seq[i].insert;
                c=(_seq[i].read == r ? '|' : 'X');
            }
            os << c;
        }
        os << "\n";
        os << "ref:    ";
        for (unsigned i(0); i<ss; ++i) os << _seq[i].ref;
        os << "\n";
        os << "insert: ";
        for (unsigned i(0); i<ss; ++i) os << _seq[i].insert;
        os << "\n";
    }

private:
    std::vector<align_position> _seq;
};
#endif


// score a contiguous matching alignment segment
//
// note that running the lnp value through as a reference creates more
// floating point stability for ambiguous alignments which have the
// same score by definition.
//
static
void
score_segment(const starling_options& /*opt*/,
              const unsigned seg_length,
              const bam_seq_base& seq,
              const uint8_t* qual,
              const unsigned read_offset,
              const bam_seq_base& ref,
              const pos_t ref_head_pos,
              double& lnp)
{

    static const double lnthird(-std::log(3.));

    for (unsigned i(0); i<seg_length; ++i)
    {
        const pos_t readi(static_cast<pos_t>(read_offset+i));
        const uint8_t sbase(seq.get_code(readi));
        if (sbase == BAM_BASE::ANY) continue;
        const uint8_t qscore(qual[readi]);
        bool is_ref(sbase == BAM_BASE::REF);
        if (! is_ref)
        {
            const pos_t refi(ref_head_pos+static_cast<pos_t>(i));
            is_ref=(sbase == ref.get_code(refi));
        }
        lnp += ( is_ref ?
                 qphred_to_ln_comp_error_prob(qscore) :
                 qphred_to_ln_error_prob(qscore)+lnthird );
    }
}



double
score_candidate_alignment(const starling_options& opt,
                          const indel_buffer& ibuff,
                          const read_segment& rseg,
                          const candidate_alignment& cal,
                          const reference_contig_segment& ref)
{
    using namespace ALIGNPATH;

#ifdef DEBUG_SCORE
    static const char GAP('-');
    align_printer ap;
#endif

    double al_lnp(0.);
    const rc_segment_bam_seq ref_bseq(ref);
    const bam_seq read_bseq(rseg.get_bam_read());
    const uint8_t* qual(rseg.qual());

    const path_t& path(cal.al.path);

#ifdef DEBUG_SCORE
    log_os << "LLAMA: path: " << path << "\n";
#endif

    unsigned read_offset(0);
    pos_t ref_head_pos(cal.al.pos);

    const std::pair<unsigned,unsigned> ends(get_match_edge_segments(path));
    const unsigned aps(path.size());
    unsigned path_index(0);
    while (path_index<aps)
    {
        const bool is_swap_start(is_segment_swap_start(path,path_index));

        unsigned n_seg(1); // number of path segments consumed
        const path_segment& ps(path[path_index]);

#ifdef DEBUG_SCORE
        log_os << "LLAMA: path_index: " << path_index << " read_offset: " << read_offset << " ref_head_pos: " << ref_head_pos << "\n";
#endif

        if       (is_swap_start)
        {
            const swap_info sinfo(path,path_index);
            n_seg=sinfo.n_seg;

            indel_key ik(ref_head_pos,INDEL::SWAP,sinfo.insert_length,sinfo.delete_length);

            // check if this is an edge swap:
            if ((path_index<ends.first) || (path_index>ends.second))
            {
                if (path_index<ends.first)
                {
                    ik=cal.leading_indel_key;
                }
                else
                {
                    ik=cal.trailing_indel_key;
                }
                assert(ik.type!=INDEL::NONE);
            }

            const indel_data* id_ptr(ibuff.get_indel_data_ptr(ik));
            if (NULL == id_ptr)
            {
                std::ostringstream oss;
                oss << "ERROR: candidate alignment does not contain expected swap indel: " << ik << "\n"
                    << "\tcandidate alignment: " << cal << "\n";
                throw blt_exception(oss.str().c_str());
            }

            const string_bam_seq insert_bseq(id_ptr->get_insert_seq());

            // if this is a leading edge-insertion we need to set
            // insert_seq_head_pos accordingly:
            //
            pos_t insert_seq_head_pos(0);
            if (path_index<ends.first)
            {
                insert_seq_head_pos=static_cast<int>(insert_bseq.size())-static_cast<int>(ps.length);
            }

            score_segment(opt,
                          sinfo.insert_length,
                          read_bseq,
                          qual,
                          read_offset,
                          insert_bseq,
                          insert_seq_head_pos,
                          al_lnp);

#ifdef DEBUG_SCORE
            for (unsigned ii(0); ii<sinfo.insert_length; ++ii)
            {
                ap.push(read_bseq.get_char(static_cast<pos_t>(read_offset+ii)),
                        GAP,
                        insert_bseq.get_char(insert_seq_head_pos+static_cast<pos_t>(ii)));
            }
            for (unsigned ii(0); ii<sinfo.delete_length; ++ii)
            {
                ap.push(GAP,
                        ref_bseq.get_char(ref_head_pos+static_cast<pos_t>(ii)),
                        GAP);
            }
#endif

        }
        else if (is_segment_align_match(ps.type))
        {
            score_segment(opt,
                          ps.length,
                          read_bseq,
                          qual,
                          read_offset,
                          ref_bseq,
                          ref_head_pos,
                          al_lnp);
#ifdef DEBUG_SCORE
            for (unsigned ii(0); ii<ps.length; ++ii)
            {
                ap.push(read_bseq.get_char(static_cast<pos_t>(read_offset+ii)),
                        ref_bseq.get_char(ref_head_pos+static_cast<pos_t>(ii)),
                        GAP);
            }
#endif

        }
        else if (ps.type==INSERT)
        {

            indel_key ik(ref_head_pos,INDEL::INSERT,ps.length);

            // check if this is an edge insertion:
            if ((path_index<ends.first) || (path_index>ends.second))
            {
                if (path_index<ends.first)
                {
                    ik=cal.leading_indel_key;
                }
                else
                {
                    ik=cal.trailing_indel_key;
                }
                assert(ik.type!=INDEL::NONE);
            }

            const indel_data* id_ptr(ibuff.get_indel_data_ptr(ik));
            if (NULL == id_ptr)
            {
                std::ostringstream oss;
                oss << "ERROR: candidate alignment does not contain expected insertion: " << ik << "\n"
                    << "\tcandidate alignment: " << cal << "\n";
                throw blt_exception(oss.str().c_str());
            }

            const string_bam_seq insert_bseq(id_ptr->get_insert_seq());

            // if this is a leading edge-insertion we need to set
            // insert_seq_head_pos accordingly:
            //
            pos_t insert_seq_head_pos(0);
            if (path_index<ends.first)
            {
                insert_seq_head_pos=static_cast<int>(insert_bseq.size())-static_cast<int>(ps.length);
            }

            score_segment(opt,
                          ps.length,
                          read_bseq,
                          qual,
                          read_offset,
                          insert_bseq,
                          insert_seq_head_pos,
                          al_lnp);

#ifdef DEBUG_SCORE
            for (unsigned ii(0); ii<ps.length; ++ii)
            {
                ap.push(read_bseq.get_char(static_cast<pos_t>(read_offset+ii)),
                        GAP,
                        insert_bseq.get_char(insert_seq_head_pos+static_cast<pos_t>(ii)));
            }
#endif

        }
        else if ((ps.type==DELETE) || (ps.type==SKIP))
        {
            // no read segment to worry about in this case
            //
#ifdef DEBUG_SCORE
            for (unsigned ii(0); ii<ps.length; ++ii)
            {
                ap.push(GAP,
                        ref_bseq.get_char(ref_head_pos+static_cast<pos_t>(ii)),
                        GAP);
            }
#endif

        }
        else if (ps.type==SOFT_CLIP)
        {
            // we rely on candidate alignment generator to suppress
            // soft-clipping so this routine does not penalizing
            // soft-clip states for now... the complication is that a
            // soft-clip alignment will always do better than its
            // unclipped equivalent. The rationale right now is that
            // if a user has soft-clipping on their input reads, they
            // want it to stay there.
            //

            // static const double lnrandom(std::log(0.25));
            // al_lnp += (ps.length*lnrandom);

        }
        else if (ps.type==HARD_CLIP)
        {
            // do nothing

        }
        else
        {
            std::ostringstream oss;
            oss << "Can't handle cigar code: " << segment_type_to_cigar_code(ps.type) << "\n";
            throw blt_exception(oss.str().c_str());
        }

        for (unsigned i(0); i<n_seg; ++i)
        {
            increment_path(path,path_index,read_offset,ref_head_pos);
        }
    }

#ifdef DEBUG_SCORE
    ap.dump(log_os);
#endif
    return al_lnp;
}

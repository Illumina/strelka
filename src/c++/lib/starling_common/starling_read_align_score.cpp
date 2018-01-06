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


/// \brief Score an insertion segment.
///
/// Note that running the \p lnp parameter through as a reference creates more
/// floating point stability compared to returning the score of the insertion segment only.
/// This stability helps us identify ambiguous alignments which have exactly the same score.
///
static
void
scoreInsertSegment(const starling_base_options& /*opt*/,
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
        const pos_t readPos(static_cast<pos_t>(read_offset+i));
        const uint8_t sbase(seq.get_code(readPos));
        if (sbase == BAM_BASE::ANY) continue;
        const uint8_t qscore(qual[readPos]);
        bool is_ref(sbase == BAM_BASE::REF);
        if (! is_ref)
        {
            const pos_t refPos(ref_head_pos+static_cast<pos_t>(i));
            is_ref=(sbase == ref.get_code(refPos));
        }
        lnp += ( is_ref ?
                 qphred_to_ln_comp_error_prob(qscore) :
                 qphred_to_ln_error_prob(qscore)+lnthird );
    }
}



/// \brief Score a match segment.
static
void
scoreMatchSegment(const unsigned seg_length,
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
        const pos_t readPos(static_cast<pos_t>(read_offset+i));
        const uint8_t sbase(seq.get_code(readPos));
        if (sbase == BAM_BASE::ANY) continue;
        const uint8_t qscore(qual[readPos]);
        bool is_ref(sbase == BAM_BASE::REF);
        if (! is_ref)
        {
            const pos_t refPos(ref_head_pos+static_cast<pos_t>(i));
            is_ref=(sbase == ref.get_code(refPos));
        }
        lnp += ( is_ref ?
                 qphred_to_ln_comp_error_prob(qscore) :
                 qphred_to_ln_error_prob(qscore)+lnthird );
    }
}



/// Convert a particular gap in the candidate alignment to the corresponding variant key to
/// lookup the variant in strelka's variant buffer.
///
static
IndelKey
getMatchingIndelKey(
    const CandidateAlignment& cal,
    const pos_t ref_head_pos,
    const unsigned delete_length,
    const unsigned insert_length,
    const std::pair<unsigned,unsigned>& ends,
    const unsigned path_index)
{
    // check if this is an edge swap:
    IndelKey indelKey;
    if ((path_index<ends.first) || (path_index>ends.second))
    {
        if (path_index<ends.first)
        {
            indelKey = cal.leading_indel_key;
        }
        else
        {
            indelKey = cal.trailing_indel_key;
        }
    }
    else
    {
        // find the indel corresponding to this point in the CIGAR alignment:
        //
        // do a dumb linear search for now...
        bool isFound(false);
        const indel_set_t& calIndels(cal.getIndels());
        for (const IndelKey& calIndel : calIndels)
        {
            if ((calIndel.pos == ref_head_pos) and
                ((calIndel.type == INDEL::INDEL) || calIndel.isMismatch()) and
                (calIndel.delete_length() == delete_length) and
                (calIndel.insert_length() == insert_length))
            {
                assert(not isFound);
                indelKey=calIndel;
                isFound=true;
            }
            else
            {
                if (calIndel.pos > ref_head_pos) break;
            }
        }
        assert(isFound);
    }

    assert(indelKey.type != INDEL::NONE);
    return indelKey;
}



/// Retrieve insertion sequence from either a breakpoint or a regular insertion.
static
const std::string&
getInsertSeq(
    const IndelKey& indelKey,
    const IndelBuffer& indelBuffer,
    const CandidateAlignment& cal)
{
    if (indelKey.is_breakpoint())
    {
        const IndelData* indelDataPtr(indelBuffer.getIndelDataPtr(indelKey));
        if (nullptr == indelDataPtr)
        {
            std::ostringstream oss;
            oss << "Candidate alignment does not contain expected breakpoint insertion: " << indelKey << "\n"
                << "\tcandidate alignment: " << cal << "\n";
            throw blt_exception(oss.str().c_str());
        }
        return indelDataPtr->getBreakpointInsertSeq();
    }
    else
    {
        return indelKey.insertSequence;
    }
}



double
scoreCandidateAlignment(
    const starling_base_options& opt,
    const IndelBuffer& indelBuffer,
    const read_segment& readSegment,
    const CandidateAlignment& cal,
    const reference_contig_segment& ref)
{
    using namespace ALIGNPATH;

#ifdef DEBUG_SCORE
    static const char GAP('-');
    align_printer ap;
#endif

    double alignmentLogProb(0.);
    const rc_segment_bam_seq ref_bseq(ref);
    const bam_seq read_bseq(readSegment.get_bam_read());
    const uint8_t* qual(readSegment.qual());

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

        IndelKey indelKey(IndelKey::noIndel());

        if       (is_swap_start || (ps.type == SEQ_MISMATCH))
        {
            unsigned deleteLength, insertLength;

            if (ps.type == SEQ_MISMATCH)
            {
                deleteLength = ps.length;
                insertLength = ps.length;
            }
            else
            {
                const swap_info sinfo(path,path_index);
                deleteLength = sinfo.delete_length;
                insertLength = sinfo.insert_length;
                n_seg=sinfo.n_seg;
            }

            indelKey = getMatchingIndelKey(cal,ref_head_pos,deleteLength,insertLength,
                                           ends,path_index);

            // a combined insert/delete event should not produce a breakpoint:
            assert(not indelKey.is_breakpoint());

            const string_bam_seq insert_bseq(indelKey.insertSequence);

            // if this is a leading edge-insertion we need to set
            // insert_seq_head_pos accordingly:
            //
            pos_t insert_seq_head_pos(0);
            if (path_index<ends.first)
            {
                insert_seq_head_pos=static_cast<int>(insert_bseq.size())-static_cast<int>(ps.length);
            }

            scoreInsertSegment(opt,
                               insertLength,
                               read_bseq,
                               qual,
                               read_offset,
                               insert_bseq,
                               insert_seq_head_pos,
                               alignmentLogProb);

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
            scoreMatchSegment(ps.length,
                              read_bseq,
                              qual,
                              read_offset,
                              ref_bseq,
                              ref_head_pos,
                              alignmentLogProb);
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
            indelKey = getMatchingIndelKey(cal,ref_head_pos,0,ps.length,
                                           ends,path_index);

            const string_bam_seq insert_bseq(getInsertSeq(indelKey,indelBuffer,cal));

            // if this is a leading edge-insertion we need to set
            // insert_seq_head_pos accordingly:
            //
            pos_t insert_seq_head_pos(0);
            if (path_index<ends.first)
            {
                insert_seq_head_pos=static_cast<int>(insert_bseq.size())-static_cast<int>(ps.length);
            }

            scoreInsertSegment(opt,
                               ps.length,
                               read_bseq,
                               qual,
                               read_offset,
                               insert_bseq,
                               insert_seq_head_pos,
                               alignmentLogProb);

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
            // no read segment to worry about in this case, but penalize non-candidate indels
            //
            if (ps.type==DELETE)
            {
                indelKey = getMatchingIndelKey(cal,ref_head_pos,ps.length,0,
                                               ends,path_index);
            }

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
            // Add a small penalty to soft-clipped segments. If this was not done, a soft-clipped alignment
            // would always have a higher likelihood than a non-clipped one, which doesn't support our intent.
            //
            // To restate: In this function we're approximately computing log of P( read | reference ) -- this
            // would be true if no candidate haplotype information was used here. Soft-clipping represents the
            // absence of an alignment between read and reference, so if the entire read were soft-clipped, then
            // P (read | reference) would be just P (read). If we assume equiprobable bases in the 'reference-vacuum'
            // case, then P(read) = 0.25^{read_length}
            //
            // The soft-clipping penalty below approximates this unaligned prior on the read for the unaligned portion
            // of the read.
            //
            static const double unalignedBasecallLogLikelihood(std::log(0.25));
            alignmentLogProb += (ps.length*unalignedBasecallLogLikelihood);
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

        // If this path segment is an indel, check if it's a candidate.
        // If it is not a candidate, add a penalty score
        // amounting to the probability that the indel occurs by chance.
        if (indelKey.type != INDEL::NONE)
        {
            const auto* indelDataPtr(indelBuffer.getIndelDataPtr(indelKey));
            assert(indelDataPtr != nullptr);
            const bool isCandidate(indelBuffer.isCandidateIndel(indelKey, *indelDataPtr));
            if (not isCandidate)
            {
                // add penalty for non-candidate indels
#if 0
                const auto& errorRates(indelDataPtr->getErrorRates());
                alignmentLogProb += errorRates.refToIndelErrorProb.getLogValue();
#else
                static const double nonCandidateIndelPenalty(std::log(1e-5));
                alignmentLogProb += nonCandidateIndelPenalty;
#endif
            }
        }

        for (unsigned i(0); i<n_seg; ++i)
        {
            increment_path(path,path_index,read_offset,ref_head_pos);
        }
    }

#ifdef DEBUG_SCORE
    ap.dump(log_os);
#endif
    return alignmentLogProb;
}

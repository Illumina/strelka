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

#include "normalizeAlignment.hh"

#include "alignment_util.hh"

#include "blt_util/align_path.hh"
#include "htsapi/align_path_bam_util.hh"

#include <cassert>



struct AlignmentInfo
{
    pos_t refPos = 0;
    pos_t readPos = 0;

    unsigned startPriorMatchSegment = 0;
    unsigned endPriorMatchSegment = 0;
    unsigned startPriorIndelSegment = 0;
    unsigned endPriorIndelSegment = 0;
    pos_t priorMatchLength = 0;
    pos_t priorDeleteLength = 0;
    pos_t priorInsertLength = 0;

    bool isChanged = false;
};



/// left-shift the current indel as detailed by the
/// contextual position values above
static
pos_t
findLeftShift(
    const bam_seq_base& refSeq,
    const bam_seq_base& readSeq,
    const AlignmentInfo& ai)
{
    pos_t shift(0);

    pos_t refRightPos(ai.refPos-1);
    pos_t readRightPos(ai.readPos-1);
    pos_t refLeftPos(refRightPos-ai.priorDeleteLength);
    pos_t readLeftPos(readRightPos-ai.priorInsertLength);

    while (shift<ai.priorMatchLength)
    {
        assert((readLeftPos-shift) >= 0);
        assert((readRightPos-shift) >= 0);
        assert((refLeftPos-shift) >= 0);
        assert((refRightPos-shift) >= 0);
        const bool isLeftMatch(readSeq.get_char(readLeftPos-shift)==refSeq.get_char(refLeftPos-shift));
        const bool isRightMatch(readSeq.get_char(readRightPos-shift)==refSeq.get_char(refRightPos-shift));
        if (isLeftMatch!=isRightMatch) break;
        shift++;
    }

    return shift;
}



static
void
leftShiftIndel(
    const bam_seq_base& refSeq,
    const bam_seq_base& readSeq,
    const unsigned currentSegment,
    alignment& al,
    AlignmentInfo& ai)
{
    using namespace ALIGNPATH;

    const pos_t shiftSize(findLeftShift(refSeq,readSeq,ai));
    if (shiftSize <= 0) return;

    ai.isChanged=true;
    assert(ai.priorMatchLength >= shiftSize);
    al.path[currentSegment].type = MATCH;
    al.path[currentSegment].length += shiftSize;
    for (unsigned matchSegIndex(ai.startPriorMatchSegment); matchSegIndex<ai.endPriorMatchSegment; ++matchSegIndex)
    {
        const unsigned newLength((matchSegIndex == ai.startPriorMatchSegment) ? (ai.priorMatchLength-shiftSize) : 0);
        al.path[matchSegIndex] = path_segment(MATCH, newLength);
    }
}



/// Indicate how much of the inserted sequence from a combined insert/delete can be matched back to the reference
/// on the left or right side.
///
/// give priority to the right side in cases that could be assigned to either, to be consistent with left-shifting
///
static
void
findInsertCollapse(
    const bam_seq_base& refSeq,
    const bam_seq_base& readSeq,
    const AlignmentInfo& ai,
    pos_t& leftCollapse,
    pos_t& rightCollapse)
{
    leftCollapse=0;
    rightCollapse=0;

    const pos_t maxCollapse(std::min(ai.priorDeleteLength,ai.priorInsertLength));

    // first try right side, then left:
    {
        const pos_t refRightPos(ai.refPos-1);
        const pos_t readRightPos(ai.readPos-1);
        while (rightCollapse < maxCollapse)
        {
            const bool isRightMatch(readSeq.get_char(readRightPos-rightCollapse)==refSeq.get_char(refRightPos-rightCollapse));
            if (! isRightMatch) break;
            rightCollapse++;
        }
    }
    {
        const pos_t refLeftPos(ai.refPos-ai.priorDeleteLength);
        const pos_t readLeftPos(ai.readPos-ai.priorInsertLength);
        while ((leftCollapse+rightCollapse) < maxCollapse)
        {
            const bool isLeftMatch(readSeq.get_char(readLeftPos+leftCollapse)==refSeq.get_char(refLeftPos+leftCollapse));
            if (! isLeftMatch) break;
            leftCollapse++;
        }
    }
}



static
void
collapseInsert(
    const bam_seq_base& refSeq,
    const bam_seq_base& readSeq,
    const unsigned currentSegment,
    alignment& al,
    AlignmentInfo& ai)
{
    using namespace ALIGNPATH;

    pos_t leftCollapse,rightCollapse;
    findInsertCollapse(refSeq,readSeq,ai,leftCollapse,rightCollapse);
    if ((leftCollapse+rightCollapse) <= 0) return;

    ai.isChanged = true;
    if (leftCollapse>0)
    {
        al.path[ai.startPriorMatchSegment].type = MATCH;
        al.path[ai.startPriorMatchSegment].length += leftCollapse;
        ai.priorMatchLength += leftCollapse;
    }
    if (rightCollapse>0)
    {
        al.path[currentSegment].type = MATCH;
        al.path[currentSegment].length += rightCollapse;
    }

    ai.priorInsertLength -= (leftCollapse+rightCollapse);
    ai.priorDeleteLength -= (leftCollapse+rightCollapse);

    bool isFirstInsert(true);
    bool isFirstDelete(true);
    for (unsigned indelSegIndex(ai.startPriorIndelSegment); indelSegIndex<ai.endPriorIndelSegment; ++indelSegIndex)
    {
        auto& ps2(al.path[indelSegIndex]);
        if (ps2.type == INSERT)
        {
            ps2.length = (isFirstInsert ? ai.priorInsertLength : 0);
            isFirstInsert = false;
        }
        else if (ps2.type == DELETE)
        {
            ps2.length = (isFirstDelete ? ai.priorDeleteLength : 0);
            isFirstDelete = false;
        }
    }
}


static
void
leftEdgeIndelCollapse(
    const bam_seq_base& refSeq,
    const bam_seq_base& readSeq,
    const unsigned currentSegment,
    alignment& al,
    AlignmentInfo& ai)
{
    using namespace ALIGNPATH;

    if (ai.priorDeleteLength > 0)
    {
        ai.isChanged = true;
        al.pos += ai.priorDeleteLength;
    }

    if (ai.priorInsertLength > 0)
    {
        pos_t rightCollapse(0);

        const pos_t refRightPos(ai.refPos-1);
        const pos_t readRightPos(ai.readPos-1);
        while ((rightCollapse < ai.priorInsertLength) && (refRightPos>0))
        {
            const bool isRightMatch(readSeq.get_char(readRightPos-rightCollapse)==refSeq.get_char(refRightPos-rightCollapse));
            if (! isRightMatch) break;
            rightCollapse++;
        }

        if (rightCollapse > 0)
        {
            ai.isChanged = true;
            al.pos -= rightCollapse;
            al.path[currentSegment].type = MATCH;
            al.path[currentSegment].length += rightCollapse;
            ai.priorInsertLength -= rightCollapse;
        }
    }

    bool isFirstInsert(true);
    for (unsigned indelSegIndex(ai.startPriorIndelSegment); indelSegIndex<ai.endPriorIndelSegment; ++indelSegIndex)
    {
        auto& ps2(al.path[indelSegIndex]);
        if (ps2.type == INSERT)
        {
            ps2.length = (isFirstInsert ? ai.priorInsertLength : 0);
            isFirstInsert = false;
        }
        else if (ps2.type == DELETE)
        {
            ps2.length = 0;
        }
    }
}



static
void
rightEdgeIndelCollapse(
    const bam_seq_base& refSeq,
    const bam_seq_base& readSeq,
    alignment& al,
    AlignmentInfo& ai)
{
    using namespace ALIGNPATH;

    if (ai.priorDeleteLength > 0)
    {
        ai.isChanged = true;
    }

    if (ai.priorInsertLength > 0)
    {
        pos_t leftCollapse(0);
        const pos_t refLeftPos(ai.refPos-ai.priorDeleteLength);
        const pos_t readLeftPos(ai.readPos-ai.priorInsertLength);
        while (leftCollapse < ai.priorInsertLength)
        {
            const bool isLeftMatch(readSeq.get_char(readLeftPos+leftCollapse)==refSeq.get_char(refLeftPos+leftCollapse));
            if (! isLeftMatch) break;
            leftCollapse++;
        }

        if (leftCollapse > 0)
        {
            ai.isChanged = true;
            al.path[ai.startPriorMatchSegment].type = MATCH;
            al.path[ai.startPriorMatchSegment].length += leftCollapse;
            ai.priorInsertLength -= leftCollapse;
        }
    }

    bool isFirstInsert(true);
    for (unsigned indelSegIndex(ai.startPriorIndelSegment); indelSegIndex<ai.endPriorIndelSegment; ++indelSegIndex)
    {
        auto& ps2(al.path[indelSegIndex]);
        if (ps2.type == INSERT)
        {
            ps2.length = (isFirstInsert ? ai.priorInsertLength : 0);
            isFirstInsert = false;
        }
        else if (ps2.type == DELETE)
        {
            ps2.length = 0;
        }
    }
}



bool
normalizeAlignment(
    const bam_seq_base& refSeq,
    const bam_seq_base& readSeq,
    alignment& al)
{
    using namespace ALIGNPATH;

    const unsigned as(al.path.size());

    //
    // first pass through is to left shift:
    //
    AlignmentInfo ai;
    {
        bool isInsideIndel(false);

        ai.refPos = al.pos;

        for (unsigned segmentIndex(0); segmentIndex<as; segmentIndex++)
        {
            auto& ps(al.path[segmentIndex]);

            if (is_segment_align_match(ps.type))
            {
                if (ai.priorMatchLength>0)
                {
                    if (isInsideIndel)
                    {
                        leftShiftIndel(refSeq,readSeq,segmentIndex,al,ai);

                        ai.priorMatchLength=0;
                    }
                }

                if (ai.priorMatchLength==0)
                {
                    ai.startPriorMatchSegment=segmentIndex;
                }
                isInsideIndel=false;
                ai.priorMatchLength += ps.length;
                ai.endPriorMatchSegment = segmentIndex+1;
            }
            else if ((ps.type == DELETE) || (ps.type == INSERT))
            {
                if (! isInsideIndel)
                {
                    isInsideIndel=true;
                    ai.priorDeleteLength=0;
                    ai.priorInsertLength=0;
                }

                if     (ps.type == DELETE)
                {
                    ai.priorDeleteLength += ps.length;
                }
                else if (ps.type == INSERT)
                {
                    ai.priorInsertLength += ps.length;
                }
            }
            else
            {
                isInsideIndel=false;
                ai.priorMatchLength=0;
            }

            if (is_segment_type_read_length(ps.type))
            {
                ai.readPos += ps.length;
            }
            if (is_segment_type_ref_length(ps.type))
            {
                ai.refPos += ps.length;
            }
        }
    }

    //
    // second pass through is to collapse internal indels:
    //
    {
        bool isInsideIndel(false);

        ai.refPos = al.pos;
        ai.readPos = 0;
        ai.priorMatchLength = 0;

        for (unsigned segmentIndex(0); segmentIndex<as; segmentIndex++)
        {
            auto& ps(al.path[segmentIndex]);
            if (ps.length == 0) continue;

            if (is_segment_align_match(ps.type))
            {
                if (ai.priorMatchLength>0)
                {
                    if (isInsideIndel)
                    {
                        if ((ai.priorDeleteLength>0) && (ai.priorInsertLength>0))
                        {
                            // collapse matching insertion sequence
                            collapseInsert(refSeq,readSeq,segmentIndex,al,ai);
                        }

                        ai.priorMatchLength=0;
                    }
                }

                if (ai.priorMatchLength==0)
                {
                    ai.startPriorMatchSegment=segmentIndex;
                }
                isInsideIndel=false;
                ai.priorMatchLength += ps.length;
                ai.endPriorMatchSegment = segmentIndex+1;
            }
            else if ((ps.type == DELETE) || (ps.type == INSERT))
            {
                if (! isInsideIndel)
                {
                    isInsideIndel=true;
                    ai.priorDeleteLength=0;
                    ai.priorInsertLength=0;
                    ai.startPriorIndelSegment=segmentIndex;
                }
                ai.endPriorIndelSegment = segmentIndex+1;

                if     (ps.type == DELETE)
                {
                    ai.priorDeleteLength += ps.length;
                }
                else if (ps.type == INSERT)
                {
                    ai.priorInsertLength += ps.length;
                }
            }
            else
            {
                isInsideIndel=false;
                ai.priorMatchLength=0;
            }

            if (is_segment_type_read_length(ps.type))
            {
                ai.readPos += ps.length;
            }
            if (is_segment_type_ref_length(ps.type))
            {
                ai.refPos += ps.length;
            }
        }
    }

    //
    // third pass is to handle edge indels
    //
    {
        bool isInsideIndel(false);

        ai.refPos = al.pos;
        ai.readPos = 0;
        ai.priorMatchLength = 0;
        ai.priorDeleteLength = 0;
        ai.priorInsertLength = 0;

        bool isFirstMatch(true);
        for (unsigned segmentIndex(0); segmentIndex<as; segmentIndex++)
        {
            auto& ps(al.path[segmentIndex]);
            if (ps.length == 0) continue;

            if (is_segment_align_match(ps.type))
            {
                if (isInsideIndel && isFirstMatch)
                {
                    leftEdgeIndelCollapse(refSeq,readSeq,segmentIndex,al,ai);
                }

                if (isFirstMatch)
                {
                    isFirstMatch=false;
                }

                if (ai.priorMatchLength>0)
                {
                    if (isInsideIndel)
                    {
                        ai.priorMatchLength=0;
                    }
                }

                if (ai.priorMatchLength==0)
                {
                    ai.startPriorMatchSegment=segmentIndex;
                }
                isInsideIndel=false;
                ai.priorMatchLength += ps.length;
                ai.endPriorMatchSegment = segmentIndex+1;
            }
            else if ((ps.type == DELETE) || (ps.type == INSERT))
            {
                if (! isInsideIndel)
                {
                    isInsideIndel=true;
                    ai.priorDeleteLength=0;
                    ai.priorInsertLength=0;
                    ai.startPriorIndelSegment=segmentIndex;
                }
                ai.endPriorIndelSegment = segmentIndex+1;

                if     (ps.type == DELETE)
                {
                    ai.priorDeleteLength += ps.length;
                }
                else if (ps.type == INSERT)
                {
                    ai.priorInsertLength += ps.length;
                }
            }
            else
            {
                if ((ps.type==SOFT_CLIP) || (ps.type==HARD_CLIP))
                {
                    if (ai.priorMatchLength && isInsideIndel)
                    {
                        if ((ai.priorDeleteLength>0) || (ai.priorInsertLength>0))
                        {
                            rightEdgeIndelCollapse(refSeq,readSeq,al,ai);
                        }
                    }

                }
                isInsideIndel=false;
                ai.priorMatchLength=0;
            }

            if (is_segment_type_read_length(ps.type))
            {
                ai.readPos += ps.length;
            }
            if (is_segment_type_ref_length(ps.type))
            {
                ai.refPos += ps.length;
            }
        }

        if (ai.priorMatchLength && isInsideIndel)
        {
            if ((ai.priorDeleteLength>0) || (ai.priorInsertLength>0))
            {
                rightEdgeIndelCollapse(refSeq,readSeq,al,ai);
            }
        }
    }

    if (ai.isChanged)
    {
        apath_cleaner(al.path);
    }
    return ai.isChanged;
}



bool
normalizeBamRecordAlignment(
    const reference_contig_segment& refSeq,
    bam_record& bamRead)
{
    alignment al;
    getAlignmentFromBamRecord(bamRead,al);

    const rc_segment_bam_seq refBamSeq(refSeq);
    const bam_seq readBamSeq(bamRead.get_bam_read());

    const bool isChanged(normalizeAlignment(refBamSeq,readBamSeq,al));

    if (isChanged)
    {
        bam1_t& br(*(bamRead.get_data()));
        br.core.pos = al.pos;
        edit_bam_cigar(al.path,br);
    }
    return isChanged;
}

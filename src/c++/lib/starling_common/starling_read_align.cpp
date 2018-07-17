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

#include "alignment_util.hh"
#include "starling_read_align.hh"
#include "starling_read_align_clipper.hh"
#include "starling_read_align_score.hh"
#include "starling_read_align_score_indels.hh"

#include "blt_util/blt_exception.hh"
#include "blt_util/log.hh"
#include "blt_util/pos_range.hh"
#include "starling_common/indel_util.hh"
#include "ActiveRegionDetector.hh"

#include <cassert>
#include <cmath>
#include <cstdlib>

#include <iostream>
#include <sstream>


//#define DEBUG_ALIGN


/// information associated with each candidate indel intersecting an alignment
struct starling_align_indel_info
{
    bool is_present = false; ///< candidate indel is present in the alignment already
    bool is_remove_only = false; ///< candidate indel can be toggled off during search but not added
    bool isInOriginalAlignment = false; ///< candidate indel is present in the original alignment
};

/// Update haplotypeConstrains based on current indel
/// haplotypeConstraints is either
/// -1: invalid
/// 0: only reference is valid
/// 1: only haplotype 1 is valid
/// 2: only haplotype 2 is valid
/// 3: both haplotypes are valid
/// \return new haplotypeConstraints
static
int
getUpdatedSampleHaplotypeConstraints(
    const int haplotypeConstraints,
    const int curIndelHaplotypeId,
    const bool isCurIndelOn,
    const bool isAnyIndelOn)
{
    // already invalid
    if (haplotypeConstraints < 0) return haplotypeConstraints;

    // curIndelHaplotypeId < 0 means that
    // any haplotype having this indel is invalid
    if ((curIndelHaplotypeId < 0) && isCurIndelOn) return -1;

    // if current indel doesn't belong to any haplotype
    // don't change haplotype constraint
    if (curIndelHaplotypeId <= 0) return haplotypeConstraints;

    // determine new haplotype constrains depending on whether current indel is on or off
    int hapConstraintsFromCurIndel = isCurIndelOn ? curIndelHaplotypeId : (3-curIndelHaplotypeId);

    // curIndelHaplotypeId is either
    // 1: only haplotype 1 is valid)
    // 2: only haplotype 2 is valid) or
    // 3: both haplotypes are valid)
    switch (hapConstraintsFromCurIndel)
    {
    case 0:
    {
        // only reference is valid
        if (isAnyIndelOn)
            return -1;
        else
            return 0;
    }
    case 1:
    case 2:
    {
        // either haplotype 1 or 2 is valid
        if ((haplotypeConstraints == 3) || (haplotypeConstraints == hapConstraintsFromCurIndel))
        {
            return hapConstraintsFromCurIndel;
        }
        if (!isAnyIndelOn)
        {
            // if no indel is on, return haplotypeConstraints=0
            // (only reference is allowed)
            return 0;
        }
        else
        {
            // this indel conflicts with the existing haplotype constraints

            return -1;
        }
    }
    case 3:
    {
        if (haplotypeConstraints > 0) return haplotypeConstraints;
        // current indel is always on here
        // invalid haplotype constrains
        return -1;
    }
    default:
        assert(false);
    }
}

/// This class holds information on phased variants in a single active region
/// and checks whether combinations of variants are valid or not
class HaplotypeStatus
{
public:
    explicit HaplotypeStatus(unsigned numSamples)
        : _haplotypeConstraints(numSamples)
    {
        for (unsigned sampleIndex(0); sampleIndex<_haplotypeConstraints.size(); ++sampleIndex)
        {
            // initially all haplotypes are valid
            _haplotypeConstraints[sampleIndex] = 3;
        }
    }

    bool updateHaplotypeStatus(
        const std::vector<int>& curIndelHaplotypeIds,
        const bool isCurIndelOn)
    {
        isAnyIndelOn = isAnyIndelOn || isCurIndelOn;

        const unsigned sampleCount(_haplotypeConstraints.size());
        bool isValid(false);
        for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
        {
            const int sampleCurIndelHaplotypeId(curIndelHaplotypeIds[sampleIndex]);
            const int sampleHaplotypeConstraints(_haplotypeConstraints[sampleIndex]);
            int updatedSampleHaplotypeConstraints =
                getUpdatedSampleHaplotypeConstraints(
                    sampleHaplotypeConstraints,
                    sampleCurIndelHaplotypeId,
                    isCurIndelOn, isAnyIndelOn);

            // this haplotype is valid if it's valid in one or more samples
            if (updatedSampleHaplotypeConstraints >= 0) isValid = true;
            _haplotypeConstraints[sampleIndex] = updatedSampleHaplotypeConstraints;
        }

        return isValid;
    }

private:
    /// haplotype constraints of all samples
    std::vector<int> _haplotypeConstraints;

    /// true if one or more indel is on
    bool isAnyIndelOn = false;
};

static
std::ostream&
operator<<(std::ostream& os, const starling_align_indel_info& ii)
{
    os << "is_present: " << ii.is_present << " is_remove_only: " << ii.is_remove_only;
    return os;
}


typedef std::map<IndelKey,starling_align_indel_info> starling_align_indel_status;

typedef std::map<ActiveRegionId,HaplotypeStatus> HaplotypeStatusMap;


static
known_pos_range
getReadAlignmentZone(
    const read_segment& rseg)
{
    const alignment& al(rseg.getInputAlignment());
    assert (! al.empty());
    return get_alignment_zone(al,rseg.read_size());
}



/// Check to see if a read segment overlaps any candidate indels (but not
/// private indels) over the maximum range suggested by its input alignment
/// as proposed by the read mapper.
///
/// If at least one overlap is discovered, then the read goes into full re-alignment.
///
/// Also run a preliminary check to make sure alignment doesn't spill outside
/// of the realignment buffer (realign_pr), if so return false. This should be
/// extremely rare.
///
static
bool
check_for_candidate_indel_overlap(
    const known_pos_range realign_buffer_range,
    const read_segment& rseg,
    const IndelBuffer& indelBuffer)
{
#ifdef DEBUG_ALIGN
    std::cerr << "BUGBUG testing read_segment for indel overlap, sample_no: " << indelBuffer.get_sample_id() << " rseg: " << rseg;
#endif

    // get liberally interpreted min, max ref coordinate bounds for the input alignment;
    const known_pos_range read_range(getReadAlignmentZone(rseg));

    // read range must be completely overlapped by the realignment buffer range
    // to continue:
    if (! realign_buffer_range.is_superset_of(read_range))
    {
        return false;
    }

#ifdef DEBUG_ALIGN
    std::cerr << "VARMIT read extends: " << read_range << "\n";
#endif

    const auto indelIterPair(indelBuffer.rangeIterator(read_range.begin_pos, read_range.end_pos));
    for (auto indelIter(indelIterPair.first); indelIter!=indelIterPair.second; ++indelIter)
    {
        const IndelKey& indelKey(indelIter->first);
#ifdef DEBUG_ALIGN
        std::cerr << "VARMIT key: " << ik;
#endif
        // check if read intersects with indel breakpoint:
        if (! is_range_intersect_indel_breakpoints(read_range,indelKey)) continue;

        const IndelData& indelData(getIndelData(indelIter));
#ifdef DEBUG_ALIGN
        std::cerr << "VARMIT intersects indel: " << ik << id;
#endif

        // check if indel qualifies as candidate indel:
        if (indelBuffer.isCandidateIndel(indelKey, indelData))
        {
#ifdef DEBUG_ALIGN
            std::cerr << "VARMIT read segment intersects at least one qualifying candidate indel.\n";
#endif
            return true;
        }
    }
#ifdef DEBUG_ALIGN
    std::cerr << "VARMIT read does not intersect indel.\n";
#endif
    return false;
}


static
void
dump_indel_status(const starling_align_indel_status& ismap,
                  std::ostream& os)
{
    for (const auto& is : ismap)
    {
        os << is.first << "status: " << is.second << "\n";
    }
}



/// is an indel either a candidate indel or in at least one of the
/// discovery alignments for this read?
///
static
bool
is_usable_indel(
    const IndelBuffer& indelBuffer,
    const IndelKey& indelKey,
    const IndelData& indelData,
    const align_id_t read_id,
    const unsigned sampleId)
{
    if (indelBuffer.isCandidateIndel(indelKey, indelData)) return true;

    const IndelSampleData& indelSampleData(indelData.getSampleData(sampleId));
    return ((indelSampleData.tier1_map_read_ids.count(read_id)>0) ||
            (indelSampleData.tier2_map_read_ids.count(read_id)>0) ||
            (indelSampleData.submap_read_ids.count(read_id)>0) ||
            (indelSampleData.noise_read_ids.count(read_id)>0));
}



/// find all indels in the indel_buffer which intersect a range (and
/// meet candidacy/usability requirements)
static
void
add_indels_in_range(
    const align_id_t read_id,
    const IndelBuffer& indelBuffer,
    const known_pos_range& pr,
    const unsigned sampleId,
    starling_align_indel_status& indel_status_map,
    std::vector<IndelKey>& indel_order)
{
    const auto indelIterPair(indelBuffer.rangeIterator(pr.begin_pos, pr.end_pos));
#ifdef DEBUG_ALIGN
    std::cerr << "VARMIT CHECKING INDELS IN RANGE: " << pr << "\n";
#endif
    for (auto indelIter(indelIterPair.first); indelIter!=indelIterPair.second; ++indelIter)
    {
        const IndelKey& indelKey(indelIter->first);
        // check if read intersects with indel and indel is usable by this read:
#ifdef DEBUG_ALIGN
        std::cerr << "VARMIT INDEL CANDIDATE " << ik;
        std::cerr << "Intersect?: " << is_range_intersect_indel_breakpoints(pr,ik) << "\n";
        std::cerr << "Usable?: " <<  is_usable_indel(indelBuffer,ik,get_indel_data(i),read_id) << "\n";
        std::cerr << "Count: " << indel_status_map.count(ik) << "\n";
#endif
        // check if the indel is not intersecting or adjacent -- if neither we don't need to
        // worry about the indel at all:
        if (! is_range_adjacent_indel_breakpoints(pr,indelKey)) continue;

        // if true, this means the indel is adjacent but not intersecting:
        const bool is_remove_only(! is_range_intersect_indel_breakpoints(pr,indelKey));

#ifdef DEBUG_ALIGN
        std::cerr << "is_remove_only " << is_remove_only << "\n";
#endif

        // if indel is already present, it may be possible to promote this indel from
        // adjacent to an intersection:
        if (indel_status_map.count(indelKey))
        {
            if ((! is_remove_only) && indel_status_map[indelKey].is_remove_only)
            {
                indel_status_map[indelKey].is_remove_only = false;
            }
        }
        else
        {
            const IndelData& indelData(getIndelData(indelIter));
            if (is_usable_indel(indelBuffer, indelKey, indelData, read_id, sampleId))
            {
                indel_status_map[indelKey].is_present = false;
                indel_status_map[indelKey].is_remove_only = is_remove_only;
                indel_order.push_back(indelKey);
            }
        }
    }
}



static
void
add_path_segment(
    ALIGNPATH::path_t& apath,
    ALIGNPATH::align_t inc_type,
    pos_t& inc_pos,
    const unsigned increment)
{
    using namespace ALIGNPATH;

    apath.push_back(path_segment(inc_type,increment));
    inc_pos += increment;
}



/// construct an alignment which includes all of the indels turned on
/// in the indel set, holding the start_position fixed to the target
/// value -- indel sets should be pre-filtered for cases where an indel
/// crosses the start pos, so this is treated as an error condition:
///
/// see unit tests
static
CandidateAlignment
make_start_pos_alignment(
    const pos_t ref_start_pos,
    const pos_t read_start_pos,
    const bool is_fwd_strand,
    const unsigned read_length,
    const indel_set_t& indels)
{
    using namespace ALIGNPATH;

    assert(read_length>0);
    assert(ref_start_pos>=0);
    assert(read_start_pos>=0);

    // if true, this read contains a leading insert,swap,softclip, etc.
    // (but not a leading deletion)
    const bool is_leading_read(read_start_pos!=0);

    CandidateAlignment cal;
    cal.al.pos=ref_start_pos;
    cal.al.is_fwd_strand=is_fwd_strand;

    pos_t ref_head_pos(ref_start_pos);
    pos_t read_head_pos(read_start_pos);

    path_t& apath(cal.al.path);

    bool isPrevMismatch(false);
    for (const IndelKey& indelKey : indels)
    {
        bool isMismatch(indelKey.isMismatch());

        // don't consider indels which can't intersect the read:
        if (indelKey.right_pos() < ref_start_pos) continue;
        if (indelKey.right_pos() == ref_start_pos)
        {
            if (isMismatch) continue;
            if (! is_leading_read) continue;
        }

        // deal with leading indel, swap or right breakpoint:
        const bool is_first_intersecting_indel(apath.empty());

        if (is_first_intersecting_indel) assert(ref_head_pos==ref_start_pos);

        if (is_leading_read && is_first_intersecting_indel)
        {
            if (indelKey.pos != ref_start_pos)
            {
                std::ostringstream oss;
                oss << "Anomalous condition for indel candidate: " << indelKey << "\n"
                    << "\tref_start_pos: " << ref_start_pos << "\n"
                    << "\tread_start_pos: " << read_start_pos << "\n"
                    << "\tref_head_pos: " << ref_head_pos << "\n"
                    << "\tread_head_pos: " << read_head_pos << "\n"
                    << "\tis_fwd_strand: " << is_fwd_strand << "\n"
                    << "\tread_length: " << read_length << "\n";

                oss << "\tfull indel set: ";
                for (const IndelKey& indelKey2 : indels)
                {
                    oss << "\t\t" << indelKey2;
                }

                throw blt_exception(oss.str().c_str());
            }

            assert((apath.empty()) && (ref_head_pos==ref_start_pos));
            assert((indelKey.insert_length() > 0) ||
                   (indelKey.type == INDEL::BP_RIGHT));

            if (indelKey.insert_length() > 0)
            {
                assert(static_cast<pos_t>(indelKey.insert_length())>=read_start_pos);
            }

            apath.push_back(path_segment(INSERT,read_start_pos));
            if (indelKey.delete_length() > 0)
            {
                add_path_segment(apath,DELETE,ref_head_pos,indelKey.delete_length());
            }
            cal.leading_indel_key=indelKey;
            isPrevMismatch = isMismatch;
            continue;
        }

        // no more leading insertion indels -- deal with regular case:
        assert((! is_first_intersecting_indel) || (read_start_pos==0));

        // note this relies on the single extra base of separation
        // required between indels during indel conflict detection:
        const bool is_edge_delete((indelKey.isPrimitiveDeletionAllele()) && (indelKey.pos == ref_start_pos));
        const int matchSegmentSize(indelKey.pos-ref_head_pos);
        const int minMatchSegmentSize((isPrevMismatch || isMismatch) ? 0 : 1);
        if ((matchSegmentSize < minMatchSegmentSize) && (!is_edge_delete))
        {
            std::ostringstream oss;
            oss << "Indel candidate: " << indelKey << " is not greater than ref_head_pos: " << ref_head_pos
                << ". ref_start_pos: " << ref_start_pos << ". Cannot resolve indel with candidate read alignment: " << cal << "\n";
            throw blt_exception(oss.str().c_str());
        }

        assert(indelKey.pos >= ref_head_pos);
        const unsigned match_segment(matchSegmentSize);

        assert(is_first_intersecting_indel || (match_segment>0) || isMismatch || isPrevMismatch);

        // remaining read segment match segment added after indel loop:
        if (((read_head_pos+match_segment) > read_length) ||
            (((read_head_pos+match_segment) == read_length) && (not indelKey.isPrimitiveDeletionAllele())))
        {
            break;
        }

        if (match_segment > 0)
        {
            apath.push_back(path_segment(MATCH,match_segment));
            ref_head_pos += match_segment;
            read_head_pos += match_segment;
        }

        if (isMismatch)
        {
            apath.push_back(path_segment(SEQ_MISMATCH, indelKey.deletionLength));
            ref_head_pos += indelKey.deletionLength;
            read_head_pos += indelKey.deletionLength;

            if (read_head_pos >= static_cast<pos_t>(read_length))
            {
                break;
            }
        }
        else if (indelKey.type == INDEL::INDEL)
        {
            if (indelKey.delete_length() > 0)
            {
                add_path_segment(apath, DELETE, ref_head_pos, indelKey.delete_length());
            }

            if (indelKey.insert_length() > 0)
            {
                const unsigned max_insert_length(read_length - read_head_pos);
                const unsigned insert_length(std::min(indelKey.insert_length(), max_insert_length));
                add_path_segment(apath, INSERT, read_head_pos, insert_length);

                const bool is_final(indelKey.insert_length() >= max_insert_length);
                if (is_final)
                {
                    cal.trailing_indel_key = indelKey;
                    break;
                }
            }
            else
            {
                if (match_segment == 0)
                {
                    cal.leading_indel_key = indelKey;
                }
                else
                {
                    const bool is_final(read_head_pos == static_cast<pos_t>(read_length));
                    if (is_final)
                    {
                        cal.trailing_indel_key = indelKey;
                    }
                }
            }
        }
        else if (indelKey.type==INDEL::BP_LEFT)
        {
            const unsigned overhang_length(read_length-read_head_pos);
            add_path_segment(apath,INSERT,read_head_pos,overhang_length);
            cal.trailing_indel_key=indelKey;
            break;
        }
        else
        {
            std::ostringstream oss;
            oss << "Unexpected indel state: " << INDEL::get_index_label(indelKey.type) << " at: " << __FILE__  << ":" << __LINE__ ;
            throw blt_exception(oss.str().c_str());
        }
        isPrevMismatch = isMismatch;
    }

    assert(read_head_pos<=static_cast<pos_t>(read_length));
    if (read_head_pos<static_cast<pos_t>(read_length))
    {
        apath.push_back(path_segment(MATCH,(read_length-read_head_pos)));
    }

    return cal;
}


/// work backwards from end_pos to get start_pos and read_start_pos
/// when the current indel set included, and then use the
/// make_start_pos_alignment routine.
///
/// see unit tests
static
void
get_end_pin_start_pos(
    const indel_set_t& indels,
    const unsigned read_length,
    const pos_t ref_end_pos,
    const pos_t read_end_pos,
    pos_t& ref_start_pos,
    pos_t& read_start_pos)
{
    assert(read_length>0);
    assert(ref_end_pos>0);
    assert(read_end_pos>0);

    ref_start_pos=ref_end_pos;
    read_start_pos=read_end_pos;

    // if true, read contains trailing insert, swap or open-ended event:
    const bool is_trailing_read(read_end_pos != static_cast<pos_t>(read_length));

    bool is_first(true);
    bool isPrevMismatch(false);
    // having trouble with normal reverse_iterator for this data
    // structure, so reversal is done by hand:
    indel_set_t::const_iterator i(indels.end()),i_begin(indels.begin());
    while (i!=i_begin)
    {
        --i;
        const IndelKey& indelKey(*i);

        const bool isMismatch(indelKey.isMismatch());
        // check that indel actually intersects the read:
        if (indelKey.pos > ref_end_pos) continue;

        if (indelKey.pos == ref_end_pos)
        {
            if (isMismatch) continue;
            if (! is_trailing_read) continue;
        }

        const bool is_trailing_indel((!isMismatch) && (indelKey.right_pos() == ref_end_pos));

        if (is_trailing_indel)   // deal with trailing-edge insert/breakpoint case first
        {
            assert((is_first) && (ref_start_pos==ref_end_pos));
            assert((indelKey.type == INDEL::INDEL) ||
                   (indelKey.type == INDEL::BP_LEFT));

            if (indelKey.type == INDEL::INDEL)
            {
                if (indelKey.insert_length() > 0)
                {
                    assert(indelKey.insert_length() >= (read_length - read_end_pos));
                }

                ref_start_pos -= indelKey.delete_length();
            }
        }
        else     // deal with normal case:
        {
            if (is_first && (read_end_pos!=static_cast<pos_t>(read_length)))
            {
                std::ostringstream oss;
                oss << "Unexpected realignment state. is_first: " << is_first
                    << " read_end_pos: " << read_end_pos
                    << " read_length: " << read_length;
                throw blt_exception(oss.str().c_str());
            }

            // note the excluding 'equals' relationship relies on the single extra base of separation
            // required between indels during indel conflict detection:
            //const bool is_edge_delete((INDEL::DELETE == indelKey.type) && (indelKey.right_pos() == ref_end_pos));

            // new indel must end at least one base below the current ref head (otherwise it would be
            // an interfering indel):
            //
            const int matchSegmentSize(static_cast<int>(ref_start_pos - indelKey.right_pos()));
            const int minMatchSegmentSize((isPrevMismatch || isMismatch) ? 0 : 1);
            if (matchSegmentSize < minMatchSegmentSize)  //&& (! is_edge_delete)) {
            {
                std::ostringstream oss;
                oss << "Unexpected indel position: indel: " << indelKey;
                oss << "\tref_start_pos: " << ref_start_pos << " ref_end_pos: " << ref_end_pos << "\n";
                throw blt_exception(oss.str().c_str());
            }

            const unsigned match_segment(std::min(matchSegmentSize,read_start_pos));

            ref_start_pos -= match_segment;
            read_start_pos -= match_segment;

            if (read_start_pos==0) return;

            if (indelKey.type == INDEL::INDEL)
            {
                ref_start_pos -= indelKey.delete_length();
                if (indelKey.insert_length() > 0)
                {
                    if (static_cast<pos_t>(indelKey.insert_length()) >= read_start_pos) return;
                    read_start_pos -= indelKey.insert_length();
                }
            }
            else if (isMismatch)
            {
                ref_start_pos -= indelKey.delete_length();
                read_start_pos -= indelKey.delete_length();

                if (read_start_pos==0) return;
            }
            else if (indelKey.type==INDEL::BP_RIGHT)
            {
                return;
            }
            else
            {
                std::ostringstream oss;
                oss << "Unexpected indel state: " << INDEL::get_index_label(indelKey.type) << " at: " << __FILE__  << ":" << __LINE__ ;
                throw blt_exception(oss.str().c_str());
            }
        }
        is_first=false;
        isPrevMismatch = isMismatch;
    }

    assert(read_start_pos >= 0);
    ref_start_pos -= read_start_pos;
    read_start_pos = 0;
}



/// to prevent incomplete search, we must put new non-present remove_only indels at the end of the list:
static
void
sort_remove_only_indels_last(const starling_align_indel_status& indel_status_map,
                             std::vector<IndelKey>& indel_order,
                             const unsigned current_depth = 0)
{
    typedef std::vector<IndelKey>::const_iterator siter;
    const siter i_begin(indel_order.begin());
    const siter i_new(i_begin+current_depth);
    const siter i_end(indel_order.end());

    std::vector<IndelKey> indel_order2;
    for (siter i(i_begin); i!=i_new; ++i)
    {
        indel_order2.push_back(*i);
    }
    for (siter i(i_new); i!=i_end; ++i)
    {
        const starling_align_indel_info& sai(indel_status_map.find(*i)->second);
        if (  (sai.is_present || (! sai.is_remove_only))) indel_order2.push_back(*i);
    }
    for (siter i(i_new); i!=i_end; ++i)
    {
        const starling_align_indel_info& sai(indel_status_map.find(*i)->second);
        if (! (sai.is_present || (! sai.is_remove_only))) indel_order2.push_back(*i);
    }
    indel_order.swap(indel_order2);
}



struct mca_warnings
{
    bool origin_skip = false;
    bool max_toggle_depth = false;
};



static
void
add_pin_exception_info(
    const char* label,
    const unsigned depth,
    const CandidateAlignment& cal,
    const CandidateAlignment& start_cal,
    const pos_t ref_start_pos,
    const pos_t read_start_pos,
    const IndelKey& cindel,
    const indel_set_t& current_indels)
{
    log_os << "\nException caught while building " << label << "-pinned alignment candidate at depth: " << depth << "\n"
           << "\tcal: " << cal
           << "\tstart_cal: " << start_cal
           << "\tref_start_pos: " << ref_start_pos << "\n"
           << "\tread_start_pos: " << read_start_pos << "\n"
           << "this_indel: " << cindel;
    for (const IndelKey& indelKey : current_indels)
    {
        log_os << "current_indels: " << indelKey;
    }
}



static
void
addKeysToCandidateAlignment(
    const starling_align_indel_status& indel_status_map,
    CandidateAlignment& cal)
{
    const auto cal_strict_pr(getStrictAlignmentRange(cal.al));
    indel_set_t calIndels;
    for (const auto& indelVal : indel_status_map)
    {
        if (not indelVal.second.is_present) continue;
        const IndelKey& indelKey(indelVal.first);
        if (not is_range_intersect_indel_breakpoints(cal_strict_pr, indelKey)) continue;
        calIndels.insert(indelKey);
    }
    if (cal.leading_indel_key.type != INDEL::NONE) calIndels.insert(cal.leading_indel_key);
    if (cal.trailing_indel_key.type != INDEL::NONE) calIndels.insert(cal.trailing_indel_key);

    cal.setIndels(calIndels);
}

/// Get haplotype IDs of current indel in all samples
static
void
getCurIndelHaplotypeIds(
    const bool isHaplotypingEnabled,
    const unsigned curSampleId,
    const IndelKey& curIndel,
    const IndelData& curIndelData,
    const bool isCurIndelInOriginalAlignment,
    std::vector<int>& curIndelHaplotypeIds)
{
    const unsigned sampleCount(curIndelHaplotypeIds.size());

    ActiveRegionId curIndelActiveRegionId(curIndelData.activeRegionId);
    if (curIndelActiveRegionId < 0) return;

    // Get haplotype IDs of current indel in all samples
    for (unsigned sampleId(0); sampleId<sampleCount; ++sampleId)
    {
        const auto& indelSampleData(curIndelData.getSampleData(sampleId));
        int curIndelHaplotypeId(indelSampleData.haplotypeId);
        if (curIndelHaplotypeId == 0)
        {
            // current indel is valid in this sample
            bool isCurIndelValidInThisSample(
                (! isHaplotypingEnabled)
                || (indelSampleData.isHaplotypingBypassed)
                || (curIndelData.isForcedOutput));
            if (!isCurIndelValidInThisSample && (sampleId == curSampleId)
                && isCurIndelInOriginalAlignment)
                isCurIndelValidInThisSample = true;

            if (curIndel.isMismatch() && (sampleId != curSampleId))
                isCurIndelValidInThisSample = false;

            if (isCurIndelValidInThisSample)
                curIndelHaplotypeId = 0;
            else
                curIndelHaplotypeId = -1;
        }
        curIndelHaplotypeIds[sampleId] = curIndelHaplotypeId;
    }
}

/// Recursively build potential alignment paths and push them into the
/// candidate alignment set:
/// \param haplotypeStatusMap map of intersecting active region and haplotyping constraints of phased variants
static
void
candidate_alignment_search(
    const starling_base_options& opt,
    const starling_base_deriv_options& dopt,
    const align_id_t read_id,
    const unsigned read_length,
    const IndelBuffer& indelBuffer,
    const unsigned sampleId,
    const known_pos_range& realign_buffer_range,
    std::set<CandidateAlignment>& cal_set,
    mca_warnings& warn,
    starling_align_indel_status indel_status_map,
    HaplotypeStatusMap haplotypeStatusMap,
    std::vector<IndelKey> indel_order,
    const unsigned depth,
    const unsigned indelToggleDepth,
    const unsigned totalToggleDepth,
    known_pos_range read_range,
    int max_read_indel_toggle,
    const CandidateAlignment& cal)
{
#ifdef DEBUG_ALIGN
    std::cerr << "VARMIT starting MCA depth: " << depth << "\n";
    std::cerr << "\twith cal: " << cal;
#endif

    // first step is to check for new indel overlaps and extend the
    // indel_status_map as necessary:
    //
    // note that we search for new indels in the range expanded from
    // the last alignment, but include an additional base from the
    // previous range so that we correctly overlap all potential new
    // indels.
    //
    bool is_new_indels(indelToggleDepth==0);
    {
        const unsigned start_ism_size(indel_status_map.size());
        const known_pos_range pr(get_soft_clip_alignment_range(cal.al));

        // check to make sure we don't realign outside of the realign buffer:
        if (! realign_buffer_range.is_superset_of(pr)) return;

        if (pr.begin_pos < read_range.begin_pos)
        {
            add_indels_in_range(read_id, indelBuffer, known_pos_range(pr.begin_pos, read_range.begin_pos + 1), sampleId,
                                indel_status_map, indel_order);
            read_range.begin_pos = pr.begin_pos;
        }
        if (pr.end_pos > read_range.end_pos)
        {
            add_indels_in_range(read_id, indelBuffer, known_pos_range(read_range.end_pos - 1, pr.end_pos), sampleId,
                                indel_status_map, indel_order);
            read_range.end_pos = pr.end_pos;
        }

        if (! is_new_indels)
        {
            is_new_indels=(start_ism_size!=indel_status_map.size());
        }

        if (is_new_indels)
        {
            sort_remove_only_indels_last(indel_status_map,indel_order,start_ism_size);
        }
    }

    // next check for recursive termination:
    if (depth == indel_order.size())
    {
        CandidateAlignment calWithKeys(cal);
        addKeysToCandidateAlignment(indel_status_map, calWithKeys);
        cal_set.insert(std::move(calWithKeys));
        return;
    }

    if (is_new_indels)
    {
        // Check for very high candidate indel density. If found, indel
        // search toggling is turned down to the minimum level which still
        // allows simple calls (distance 1 from input alignment). The intention
        // is to allow basic indel calling to proceed in practical time
        // through regions with very high numbers of candidate indels.
        //
        // TODO: Note the expression for indel density doesn't account for
        // a possibly large number of indels intersected by a large
        // deletion in a read. This works well enough for now, but if the
        // max indel size is ever run around the order of 10k or more this
        // might start to spuriously engage the filter.
        //
        const double max_indels(read_length*opt.max_candidate_indel_density);
        if (indel_status_map.size()>max_indels)
        {
            max_read_indel_toggle=1;
        }
        else
        {
            max_read_indel_toggle=opt.max_read_indel_toggle;
        }

        // a new stronger complexity limit on search based on total candidate indels crossing the read:
        //
        {
            const int max_toggle(dopt.sal.get_max_toggle(indel_status_map.size()));
            max_read_indel_toggle=std::min(max_read_indel_toggle,max_toggle);
        }
    }

    // check whether toggling the input alignment already exceeds the maximum
    // number of toggles made to the exemplar alignment (this is
    // here to prevent a combinatorial blowup)
    //
    if (static_cast<int>(indelToggleDepth)>max_read_indel_toggle)
    {
        warn.max_toggle_depth=true;
        return;
    }

    // each recursive step invokes (up to) 3 paths:
    //  1) is the current state of the active indel
    //  2) is the alternate state of the active indel with the alignment's start position pinned
    //  3) is the alternate state of the active indel with the alignment's end position pinned
    //
    // toggle will not be invoked for unusable indels (but those are
    // already filtered out of the list)
    //
    // ??? -- no longer true??? -- toggle will not be invoked if they
    // lead to indel conflicts with a previous indel
    //
    // start or end position may be skipped if a deletion spans one or
    // both of these points
    //
    // edge-indels can only be pinned on one side
    //

    const IndelKey& curIndel(indel_order[depth]);

    bool isCurIndelConflicting(false);
    bool containsIndelNotDiscoveredFromReads(false);
    for (unsigned i(0); i<depth; ++i)
    {
        const IndelKey& indelKey(indel_order[i]);

        if (!(indel_status_map[indelKey].is_present)) continue;
        if (is_indel_conflict(indelKey, curIndel)) isCurIndelConflicting = true;
        if ((! containsIndelNotDiscoveredFromReads) &&
            (indelBuffer.getIndelDataPtr(indelKey)->status.notDiscoveredFromReads))
        {
            containsIndelNotDiscoveredFromReads = true;
        }
    }

    const unsigned sampleCount(opt.getSampleCount());

    const bool isCurIndelOn(indel_status_map[curIndel].is_present);
    const IndelData& curIndelData(*indelBuffer.getIndelDataPtr(curIndel));

    const ActiveRegionId curIndelActiveRegionId(curIndelData.activeRegionId);
    const bool isCurIndelInActiveRegion(curIndelActiveRegionId >= 0);
    if (isCurIndelInActiveRegion)
    {
        auto it = haplotypeStatusMap.find(curIndelActiveRegionId);
        if (it == haplotypeStatusMap.end())
        {
            haplotypeStatusMap.insert(std::make_pair(curIndelActiveRegionId, HaplotypeStatus(sampleCount)));
        }
    }

    // true if current indel is an external indel not discovered from reads
    const bool isCurIndelNotDiscoveredFromReads(curIndelData.status.notDiscoveredFromReads);

    // Get haplotype IDs of current indel in all samples
    std::vector<int> curIndelHaplotypeIds(sampleCount);
    getCurIndelHaplotypeIds(
        opt.isHaplotypingEnabled,
        sampleId,
        curIndel,
        curIndelData,
        indel_status_map[curIndel].isInOriginalAlignment,
        curIndelHaplotypeIds);

    // alignment 1) --> unchanged case:
    try
    {
        bool isNextCandidateAlignmentValid(true);

        // 1. Check haplotype constraints
        HaplotypeStatusMap newHaplotypeStatusMap(haplotypeStatusMap);
        if ((! isCurIndelConflicting) && isCurIndelInActiveRegion)
        {
            isNextCandidateAlignmentValid = newHaplotypeStatusMap.at(curIndelActiveRegionId)
                                            .updateHaplotypeStatus(curIndelHaplotypeIds, isCurIndelOn);
        }
        else
        {
            // current indel is conflicting or there's no haplotype info
            isNextCandidateAlignmentValid = (! curIndel.isMismatch()) || (! isCurIndelOn);
        }

        // 2. Check whether we will be including
        // more than one indels without enough read support
        if (isCurIndelOn && containsIndelNotDiscoveredFromReads &&
            isCurIndelNotDiscoveredFromReads)
        {
            // it's prohibited to include more than one indels without enough read support
            isNextCandidateAlignmentValid = false;
        }

        if ((! isNextCandidateAlignmentValid) && (totalToggleDepth == 0))
        {
            // even if the above conditions are not met,
            // if there was no indel toggle,
            // allow the alignment search to proceed
            isNextCandidateAlignmentValid = true;
        }

        if (isNextCandidateAlignmentValid)
        {
            candidate_alignment_search(opt, dopt, read_id, read_length, indelBuffer,
                                       sampleId, realign_buffer_range,
                                       cal_set, warn,
                                       indel_status_map, newHaplotypeStatusMap,
                                       indel_order, depth + 1, indelToggleDepth, totalToggleDepth,
                                       read_range, max_read_indel_toggle, cal);
        }
    }
    catch (...)
    {
        log_os << "\nException caught while building default alignment candidate at depth: " << depth << "\n"
               << "\tcal: " << cal
               << "this_indel: " << curIndel;
        throw;
    }

    bool isNextCandidateAlignmentValid(true);
    HaplotypeStatusMap newHaplotypeStatusMap(haplotypeStatusMap);
    if (!isCurIndelConflicting && isCurIndelInActiveRegion)
    {
        isNextCandidateAlignmentValid = newHaplotypeStatusMap.at(curIndelActiveRegionId).updateHaplotypeStatus(curIndelHaplotypeIds, !isCurIndelOn);
    }
    else
    {
        isNextCandidateAlignmentValid = !curIndel.isMismatch() || isCurIndelOn;
    }

    // Check whether we will be including
    // more than one indels not discovered from reads
    if (!isCurIndelOn &&
        containsIndelNotDiscoveredFromReads &&
        isCurIndelNotDiscoveredFromReads)
    {
        // it's prohibited to include more than one indels without enough read support
        isNextCandidateAlignmentValid = false;
    }

    if (! isNextCandidateAlignmentValid) return;

    if (! isCurIndelOn)
    {
        // check whether this is a remove only indel:
        if (indel_status_map[curIndel].is_remove_only) return;

        // check whether this indel would interfere with an indel that's
        // already been toggled on:
        //
        if (isCurIndelConflicting) return;
    }

    // Mismatches discovered in AR doesn't increase toggle depth,
    // because #alignments is constrained by phasing info
    unsigned indelToggleIncrement(curIndel.isMismatch() ? 0 : 1);

    // check whether toggling this indel would exceed the maximum
    // number of toggles made to the exemplar alignment (this is
    // here to prevent a combinatorial blowup)
    //
    if (static_cast<int>(indelToggleDepth+indelToggleIncrement)>max_read_indel_toggle)
    {
        warn.max_toggle_depth=true;
        return;
    }

    // changed cases:
    indel_status_map[curIndel].is_present=(! isCurIndelOn);

    // extract only those indels that are present in the next
    // alignment:
    //
    indel_set_t current_indels;
    for (const auto& is : indel_status_map)
    {
        if (is.second.is_present) current_indels.insert(is.first);
    }

    // a pin on either end of the alignment is not possible/sensible
    // if:
    //
    // A) a deletion is being added which spans the pin site
    // B) an edge insertion/breakpoint is being removed from the pinned side

    // alignment 2) -- insert or delete indel and pin the start position
    //
    // test for conditions where the start pin is not possible:
    //
    {
        const pos_t ref_start_pos(cal.al.pos);

        bool is_start_pin_valid(true);
        if (! curIndel.isMismatch())
        {
            const bool is_start_pos_delete_span(curIndel.open_pos_range().is_pos_intersect(ref_start_pos));
            const bool is_start_pos_indel_span(isCurIndelOn && (curIndel == cal.leading_indel_key));
            is_start_pin_valid = (! (is_start_pos_delete_span || is_start_pos_indel_span));
        }

#ifdef DEBUG_ALIGN
        std::cerr << "VARMIT toggling MCA depth: " << depth << "\n";
        std::cerr << "VARMIT current indel: " << cindel;
        std::cerr << "VARMIT current indel on?: " << is_cindel_on << "\n";
        std::cerr << "VARMIT start-pin valid?: " << is_start_pin_valid << "\n";
#endif

        if (is_start_pin_valid)
        {
            const pos_t read_start_pos(unalignedPrefixSize(cal.al.path));
            CandidateAlignment start_cal;
            try
            {
                start_cal = make_start_pos_alignment(ref_start_pos,
                                                     read_start_pos,
                                                     cal.al.is_fwd_strand,
                                                     read_length,
                                                     current_indels);

                candidate_alignment_search(opt, dopt, read_id, read_length, indelBuffer,
                                           sampleId,
                                           realign_buffer_range, cal_set,
                                           warn, indel_status_map, newHaplotypeStatusMap,
                                           indel_order, depth + 1, indelToggleDepth + indelToggleIncrement, totalToggleDepth + 1,
                                           read_range, max_read_indel_toggle,
                                           start_cal);
            }
            catch (...)
            {
                add_pin_exception_info("start",depth,cal,start_cal,ref_start_pos,read_start_pos,curIndel,current_indels);
                throw;
            }
        }
    }

    // check whether this is a mismatch or an equal-length swap,
    // in which case alignment 3 is unnecessary:
    if (curIndel.isMismatch()) return;
    if ((curIndel.type==INDEL::INDEL) && (curIndel.delete_length()==curIndel.insert_length())) return;

    // alignment 3) -- insert or delete indel and pin the end position
    //
    // test for conditions where end-pin is not possible:
    //
    {
        const pos_t ref_end_pos(cal.al.pos+apath_ref_length(cal.al.path));

        // end pin is not possible when
        // (1) an indel deletes through the end-pin position
        // (2) we try to remove a trailing indel [TODO seems like same rule should be in place for adding a trailing indel]
        const bool is_end_pos_delete_span(curIndel.open_pos_range().is_pos_intersect(ref_end_pos-1));
        const bool is_end_pos_indel_span(isCurIndelOn && (curIndel == cal.trailing_indel_key));
        const bool is_end_pin_valid(! (is_end_pos_delete_span || is_end_pos_indel_span));

#ifdef DEBUG_ALIGN
        std::cerr << "VARMIT toggling MCA depth: " << depth << "\n";
        std::cerr << "VARMIT current indel: " << cindel;
        std::cerr << "VARMIT current indel on?: " << is_cindel_on << "\n";
        std::cerr << "VARMIT end-pin valid?: " << is_end_pin_valid << "\n";
#endif

        if (is_end_pin_valid)
        {
            // work backwards from end_pos to get start_pos and
            // read_start_pos when the current indel set included,
            // and then used the make_start_pos_alignment routine.
            const pos_t read_end_pos(read_length-+unalignedSuffixSize(cal.al.path));
            pos_t ref_start_pos(0);
            pos_t read_start_pos(0);
            get_end_pin_start_pos(current_indels,read_length,
                                  ref_end_pos,read_end_pos,
                                  ref_start_pos,read_start_pos);

            // guard against low-frequency circular chromosome event:
            if (ref_start_pos<0)
            {
                warn.origin_skip=true;
            }
            else
            {
                CandidateAlignment start_cal;
                try
                {
                    start_cal = make_start_pos_alignment(ref_start_pos,
                                                         read_start_pos,
                                                         cal.al.is_fwd_strand,
                                                         read_length,
                                                         current_indels);

                    candidate_alignment_search(opt, dopt, read_id, read_length, indelBuffer, sampleId,
                                               realign_buffer_range, cal_set,
                                               warn, indel_status_map, newHaplotypeStatusMap,
                                               indel_order, depth + 1, indelToggleDepth + indelToggleIncrement, totalToggleDepth + 1,
                                               read_range, max_read_indel_toggle, start_cal);
                }
                catch (...)
                {
                    add_pin_exception_info("end",depth,cal,start_cal,ref_start_pos,read_start_pos,curIndel,current_indels);
                    log_os << "ref_end_pos: " << ref_end_pos << "\n"
                           << "read_end_pos: " << read_end_pos << "\n";
                    throw;
                }
            }
        }
    }
}


/// Summary stats of an alignment
struct extra_path_info
{
    unsigned indelCount = 0;
    unsigned totalDeletionSize = 0;
    unsigned totalInsertionSize = 0;

    /// Sum of CIGAR segment start positions - this is a meaningless statistic used to tie break otherwise
    /// identical alignments
    unsigned sumSegmentPos = 0;
};



static
extra_path_info
getExtraPathInfo(const ALIGNPATH::path_t& p)
{
    using namespace ALIGNPATH;

    unsigned read_pos(0);

    extra_path_info epi;
    for (const path_segment& ps : p)
    {
        if (! is_segment_align_match(ps.type)) epi.indelCount++;
        if (ps.type == DELETE)
        {
            epi.totalDeletionSize += ps.length;
            epi.sumSegmentPos += read_pos;
        }
        if (ps.type == INSERT)
        {
            epi.totalInsertionSize += ps.length;
            epi.sumSegmentPos += read_pos;
        }

        if (is_segment_type_read_length(ps.type)) read_pos += ps.length;
    }

    return epi;
}



static
unsigned int
getCandidateIndelCount(
    const IndelBuffer& indelBuffer,
    const CandidateAlignment& cal)
{
    unsigned val(0);
    for (const IndelKey& indelKey : cal.getIndels())
    {
        if (indelBuffer.isCandidateIndel(indelKey)) val++;
    }
    return val;
}



/// \brief Determine which of two alignments is the best one to use as the 'representative' one for SNV calling
///
/// It is assumed that the caller has already determined some high level of alignment score similarity
///
/// The function will select the alignment with the following properties, in order until a tie is broken:
/// 1. fewest indels
/// 2. fewest non-candidate indels
/// 2. smallest total insertion length
/// 3. smallest total deletion length
///
static
bool
isFirstCandidateAlignmentPreferred(
    const IndelBuffer& indelBuffer,
    const CandidateAlignment& c1,
    const CandidateAlignment& c2)
{
    const extra_path_info epi1(getExtraPathInfo(c1.al.path));
    const extra_path_info epi2(getExtraPathInfo(c2.al.path));

    if (epi2.indelCount < epi1.indelCount) return false;
    if (epi2.indelCount > epi1.indelCount) return true;

    const unsigned cic1(getCandidateIndelCount(indelBuffer, c1));
    const unsigned cic2(getCandidateIndelCount(indelBuffer, c2));
    if (cic2 > cic1) return false;
    if (cic2 < cic1) return true;

    if (epi2.totalInsertionSize < epi1.totalInsertionSize) return false;
    if (epi2.totalInsertionSize > epi1.totalInsertionSize) return true;

    if (epi2.totalDeletionSize < epi1.totalDeletionSize) return false;
    if (epi2.totalDeletionSize > epi1.totalDeletionSize) return true;

    return (epi2.sumSegmentPos >= epi1.sumSegmentPos);
}



#if 0
// make sure the cal pool contains at least one candidate indel:
static
bool
is_cal_pool_contains_candidate(const starling_base_options& client_opt,
                               const depth_buffer& db,
                               const indel_buffer& ibuff,
                               const cal_pool_t& max_cal_pool)
{

    const unsigned n_cal(max_cal_pool.size());
    for (unsigned i(0); i<n_cal; ++i)
    {
        if (0 < get_candidate_indel_count(client_opt,db,ibuff,*(max_cal_pool[i]))) return true;
    }
    return false;
}
#endif


/// \brief Record 'best' re-alignment for the purpose of defining a single realigned read
///
/// Note that 'best' is not always the most probable alignment. This routine will optionally
/// add soft-clipping if the top candidate alignment pool suggests ambiguous edge regions.
///
/// \param[in] topAlignmentPtrs pointers to high scoring candidate alignments
/// \param[in] bestAlignmentPtr pointer to the 'best' alignment as determined by the caller;
///          the 'best' alignment must be present in the top alignment set
/// \param[out] realignment the final realignment output by this function
static
void
finishRealignment(
    const cal_pool_t& topAlignmentPtrs,
    const CandidateAlignment* bestAlignmentPtr,
    alignment& realignment)
{
    if (topAlignmentPtrs.size() > 1)
    {
        // soft-clip off any ambiguous regions from the alignment:
        // NOTE this can result in an empty alignment!!!
        //
        const unsigned topAlignmentCount(topAlignmentPtrs.size());
        unsigned bestAlignmnetIndex(topAlignmentCount);
        for (unsigned alignmentIndex(0); alignmentIndex<topAlignmentCount; ++alignmentIndex)
        {
            if (topAlignmentPtrs[alignmentIndex]==bestAlignmentPtr)
            {
                bestAlignmnetIndex=alignmentIndex;
                break;
            }
        }
        assert(bestAlignmnetIndex != topAlignmentCount);
        getClippedAlignmentFromTopAlignmentPool(topAlignmentPtrs, bestAlignmnetIndex, realignment);
        if (realignment.empty())
        {
            realignment = bestAlignmentPtr->al;
#ifdef DEBUG_ALIGN
            log_os << "VARMIT clipping failed -- revert to: " << realignment;
        }
        else
        {
            log_os << "VARMIT clipped alignment: " << realignment;
#endif
        }
    }
    else
    {
        realignment = bestAlignmentPtr->al;
    }
}



static
bool
is_alignment_spanned_by_range(const known_pos_range pr,
                              const alignment& al)
{
    return pr.is_superset_of(getStrictAlignmentRange(al));
}



static
void
bam_seq_to_str(
    const bam_seq_base& bs,
    const unsigned start,
    const unsigned end,
    std::string& s)
{
    s.clear();
    for (unsigned i(start); i<end; ++i) s.push_back(bs.get_char(i));
}



/// Initialize a candidate alignment from a standard alignment.
///
static
void
getCandidateAlignment(
    const alignment& al,
    const read_segment& rseg,
    CandidateAlignment& cal)
{
    using namespace ALIGNPATH;

    cal.al=al;

    pos_t read_pos(0);
    pos_t ref_pos(al.pos);

    const std::pair<unsigned,unsigned> ends(get_match_edge_segments(al.path));
    const unsigned as(al.path.size());
    for (unsigned i(0); i<as; ++i)
    {
        const path_segment& ps(al.path[i]);
        if ((INSERT == ps.type) || (DELETE == ps.type))
        {
            IndelKey indelKey(ref_pos,INDEL::INDEL);
            if (INSERT == ps.type)
            {
                bam_seq_to_str(rseg.get_bam_read(),read_pos,read_pos+ps.length,indelKey.insertSequence);
            }
            else
            {
                indelKey.deletionLength = ps.length;
            }

            if     (i<ends.first)
            {
                cal.leading_indel_key = indelKey;
            }
            else if (i>ends.second)
            {
                cal.trailing_indel_key = indelKey;
            }
        }
        if (is_segment_type_read_length(ps.type)) read_pos += ps.length;
        if (is_segment_type_ref_length(ps.type)) ref_pos += ps.length;
    }
}



/// \brief Score the candidate alignments and identify/create alignments with privileged roles in downstream calling
///
/// The three operations performed here are:
/// (1) Score the candidate alignments
/// (2) Identify the highest scoring alignment
/// (3) Identify the 'representative' alignment to use for downstream SNV calling. This is chosen based on
///     a combination of high score and low complexity.
///
static
void
scoreCandidateAlignments(
    const starling_base_options& opt,
    const reference_contig_segment& ref,
    read_segment& readSegment,
    IndelBuffer& indelBuffer,
    const std::set<CandidateAlignment>& candAlignments,
    std::vector<double>& candAlignmentScores,
    double& maxCandAlignmentScore,
    const CandidateAlignment*& maxCandAlignmentPtr,
    const bool isTestSoftClippedInputAligned,
    const alignment& softClippedInputAlignment)
{
    // The smooth optimum alignment and alignment pool are
    // used for realignment, whereas the strict max_path path
    // alignment is reported back to the indel_scoring routines.
    //
    // The smooth set may be different from the max pool set for two
    // reasons: (1) the path score is within smooth range of the
    // optimum (2) the smooth set is restrained by one or both edges of
    // an alignment being 'pinned' for exon alignments.
    //

    // pins are used to prevent exon/intron boundaries from being moved
    // during exon realignment:
    //
    const std::pair<bool,bool> edge_pin(readSegment.get_segment_edge_pin());
    const bool is_pinned(edge_pin.first || edge_pin.second);

    const auto cal_set_begin(candAlignments.cbegin()), cal_set_end(candAlignments.cend());
    for (auto cal_iter(cal_set_begin); cal_iter!=cal_set_end; ++cal_iter)
    {
        const CandidateAlignment& ical(*cal_iter);
        const double path_lnp(
            scoreCandidateAlignment(opt, indelBuffer, readSegment, ical, ref));

        candAlignmentScores.push_back(path_lnp);

#ifdef DEBUG_ALIGN
        std::cerr << "VARMIT CANDIDATE ALIGNMENT " << ical;
        std::cerr << "score: " << path_lnp << "\n";
#endif

        // check that this is not the first iteration, then determine if this
        // iteration's candidate alignment will be the new max value:
        //
        if (nullptr != maxCandAlignmentPtr)
        {
            if (path_lnp<maxCandAlignmentScore) continue;

            // TODO -- cleaner test of float equivalence (the
            // present test should be legit given the way the
            // score is calculated, but it's still not preferred)
            //
            if ((path_lnp<=maxCandAlignmentScore) &&
                isFirstCandidateAlignmentPreferred(indelBuffer, *maxCandAlignmentPtr, ical)) continue;
        }
        maxCandAlignmentScore=path_lnp;
        maxCandAlignmentPtr=&ical;
    }

    // if there's a pin on this segment, then find out which
    // alignments are allowed and what the max path_lnp of this
    // subset is:
    //
    double max_allowed_path_lnp(maxCandAlignmentScore);
    std::vector<bool> is_cal_allowed;
    if (is_pinned)
    {
        const CandidateAlignment* max_allowed_cal_ptr(nullptr);
        is_cal_allowed.resize(candAlignments.size(),true);
        const known_pos_range startingAlignmentRange(getStrictAlignmentRange(readSegment.getInputAlignment()));

        unsigned cal_index(0);
        for (auto cal_iter(cal_set_begin); cal_iter!=cal_set_end; ++cal_iter,++cal_index)
        {
            const known_pos_range candidateAlignmentRange(getStrictAlignmentRange(cal_iter->al));
            if       (edge_pin.first && (candidateAlignmentRange.begin_pos != startingAlignmentRange.begin_pos))
            {
                is_cal_allowed[cal_index] = false;
            }
            else if (edge_pin.second && (candidateAlignmentRange.end_pos != startingAlignmentRange.end_pos))
            {
                is_cal_allowed[cal_index] = false;
            }
            if (! is_cal_allowed[cal_index]) continue;

            if (nullptr!=max_allowed_cal_ptr)
            {
                if (candAlignmentScores[cal_index]<max_allowed_path_lnp) continue;
            }
            max_allowed_path_lnp=candAlignmentScores[cal_index];
            max_allowed_cal_ptr=&(*cal_iter);
        }

        if (nullptr == max_allowed_cal_ptr)
        {
            std::ostringstream oss;
            oss << "Reached anomalous state during search for most likely exon alignment.\n";
            oss << "\tread_segment: " << readSegment << "\n";
            oss << "\tCandidate alignments:\n";
            for (auto cal_iter(cal_set_begin); cal_iter!=cal_set_end; ++cal_iter)
            {
                oss << *cal_iter << "\n";
            }
            throw blt_exception(oss.str().c_str());
        }
    }

    // go back through the the path_lnp values and produce a pool
    // that:
    //
    // (1) possibly allows sub maximal values (if
    // is_smoothed_alignments is set)
    //
    // (2) accounts for end pins (used on exon splice sites for
    // instance)
    //
    double allowed_lnp_range(0);
    if (opt.is_smoothed_alignments)
    {
        allowed_lnp_range=opt.smoothed_lnp_range;
    }

    double smooth_path_lnp(0);
    const CandidateAlignment* smooth_cal_ptr(nullptr);
    cal_pool_t smooth_cal_pool; // set of alignments within smooth-range of max path score

    unsigned cal_index(0);
    for (auto cal_iter(cal_set_begin); cal_iter!=cal_set_end; ++cal_iter,++cal_index)
    {
        if ((candAlignmentScores[cal_index]+allowed_lnp_range) < max_allowed_path_lnp) continue;
        if (is_pinned && (! is_cal_allowed[cal_index])) continue;
        const CandidateAlignment& ical(*cal_iter);
        smooth_cal_pool.push_back(&ical);
        if ((nullptr==smooth_cal_ptr) ||
            (!isFirstCandidateAlignmentPreferred(indelBuffer, *smooth_cal_ptr, ical)))
        {
            smooth_path_lnp=candAlignmentScores[cal_index];
            smooth_cal_ptr=&ical;
        }
    }

    assert(nullptr != smooth_cal_ptr);

#ifdef DEBUG_ALIGN
    std::cerr << "BUBBY: key,max_path_lnp,max_path: " << rseg.key() << " " << maxCandAlignmentScore << " max_cal: " << *maxCandAlignmentPtr;

    if (smooth_cal_pool.size() > 1)
    {
        const unsigned n_cal(smooth_cal_pool.size());
        std::cerr << "BUBBY: " << n_cal << " final alignment pool:\n";
        for (unsigned i(0); i<n_cal; ++i)
        {
            std::cerr << "BUBBY: alignment " << i << "\n" << *(smooth_cal_pool[i]);
            const known_pos_range ipr(get_strict_alignment_range(smooth_cal_pool[i]->al));
            for (unsigned j(i+1); j<n_cal; ++j)
            {
                const known_pos_range jpr(get_strict_alignment_range(smooth_cal_pool[j]->al));
                if (ipr.begin_pos==jpr.begin_pos && ipr.end_pos==jpr.end_pos) std::cerr << "COWSLIP\n";
            }
        }
    }

    if (opt.is_smoothed_alignments)
    {
        std::cerr << "BUBBY: smooth_path_lnp,smooth_path: " << smooth_path_lnp << " smooth_cal: " << *smooth_cal_ptr;
    }
#endif

    //
    // record the read realignment which will be used for pileup and snp-calling:
    //

    // this option allows the original read alignment (with soft-clipping) to be used for the pileup
    // if a better alignment was not found:
    if (isTestSoftClippedInputAligned)
    {
        // convert alignment into CandidateAlignment type:
        CandidateAlignment softClippedCandidateAlignment;
        {
            getCandidateAlignment(softClippedInputAlignment, readSegment, softClippedCandidateAlignment);
            const bool includeMismatches(false);
            indel_set_t candidateAlignmentIndels;
            getAlignmentIndels(softClippedCandidateAlignment, ref, readSegment,
                               opt.maxIndelSize, includeMismatches, candidateAlignmentIndels);
            softClippedCandidateAlignment.setIndels(candidateAlignmentIndels);
        }

        // score candidate alignment:
        const double path_lnp(scoreCandidateAlignment(opt, indelBuffer, readSegment,
                                                      softClippedCandidateAlignment, ref));

        if (path_lnp >= smooth_path_lnp)
        {
            if (not (softClippedInputAlignment == readSegment.getInputAlignment()))
            {
                readSegment.is_realigned = true;
                readSegment.realignment = softClippedInputAlignment;
            }
            return;
        }
    }

    readSegment.is_realigned = true;
    finishRealignment(smooth_cal_pool, smooth_cal_ptr, readSegment.realignment);
}



/// \brief Find the most likely alignment and most likely alignment for
/// each indel state for every indel in indel_status_map
///
static
void
scoreCandidateAlignmentsAndIndels(
    const starling_base_options& opt,
    const starling_base_deriv_options& dopt,
    const starling_sample_options& sample_opt,
    const reference_contig_segment& ref,
    read_segment& rseg,
    IndelBuffer& indelBuffer,
    const unsigned sampleId,
    std::set<CandidateAlignment>& candAlignments,
    const bool is_incomplete_search,
    const bool isTestSoftClippedInputAligned,
    const alignment& softClippedInputAlignment)
{
    assert(! candAlignments.empty());

    // (1) score all alignments in cal_set and find the max scoring
    // alignment path.
    //
    // NOTE: In (1) we're essentially picking a "best" haplotype by
    // searching for the best alignment with no-penalty candidate indels applied. If we also had candidate snps,
    // it would be a true candidate haplotype search.
    //
    // (1b) possibly also find a "best" scoring path for re-alignment
    // which is not (quite) the max
    //
    std::vector<double> candAlignmentScores;
    double maxCandAlignmentScore(0);
    const CandidateAlignment* maxCandAlignmentPtr(nullptr);

    try
    {
        scoreCandidateAlignments(opt, ref, rseg, indelBuffer, candAlignments,
                                 candAlignmentScores, maxCandAlignmentScore, maxCandAlignmentPtr,
                                 isTestSoftClippedInputAligned, softClippedInputAlignment);
    }
    catch (...)
    {
        log_os << "Exception caught while scoring candidate alignments for read segment: " << rseg;
        throw;
    }

    // Realignment for snp-calling and visualization is complete here,
    // remaining task is to evaluate alternate alignments as required
    // by the indel calling model.
    //

    // For snps and visualization we realign all reads that aren't filtered out entirely,
    // but for the score_indels step, we only consider reads meeting the minimum filtration
    // criteria setup for indels.
    //
    if (! rseg.is_tier1or2_mapping()) return;

    try
    {
        score_indels(opt, dopt, sample_opt, rseg, indelBuffer, sampleId, candAlignments, is_incomplete_search,
                     candAlignmentScores, maxCandAlignmentScore, maxCandAlignmentPtr);
    }
    catch (...)
    {
        log_os << "Exception caught while scoring indels for read segment: " << rseg;
        throw;
    }
}



static
void
getCandidateAlignments(
    const starling_base_options& opt,
    const starling_base_deriv_options& dopt,
    const reference_contig_segment& ref,
    const read_segment& rseg,
    const IndelBuffer& indelBuffer,
    const unsigned sampleId,
    const alignment& inputAlignment,
    const known_pos_range realign_buffer_range,
    mca_warnings& warn,
    std::set<CandidateAlignment>& cal_set)
{
    const unsigned read_length(rseg.read_size());

    starling_align_indel_status indel_status_map;
    HaplotypeStatusMap haplotypeStatusMap;
    std::vector<IndelKey> indel_order;

    CandidateAlignment cal;
    getCandidateAlignment(inputAlignment, rseg, cal);

#ifdef DEBUG_ALIGN
    std::cerr << "VARMIT starting search from input alignment: " << cal;
#endif

    // Get indel set and indel order for the input alignment:
    const known_pos_range exemplar_pr(get_soft_clip_alignment_range(cal.al));
    add_indels_in_range(rseg.getReadIndex(), indelBuffer, exemplar_pr, sampleId, indel_status_map, indel_order);

#ifdef DEBUG_ALIGN
    std::cerr << "VARMIT exemplar alignment range: " << exemplar_pr << "\n";
#endif

    // Mark the indels which are already included in the input
    // alignment.
    //
    {
        indel_set_t candidateAlignmentIndels;
        const bool includeMismatches(true);
        // get alignment indels including mismatches
        getAlignmentIndels(cal, ref, rseg, opt.maxIndelSize, includeMismatches, candidateAlignmentIndels);

        indel_set_t validIndels;
        bool isRecomputeRequired(false);
        for (const IndelKey& indelKey : candidateAlignmentIndels)
        {
            if (indel_status_map.find(indelKey)==indel_status_map.end())
            {
                if (indelKey.isMismatch()) continue;

                std::ostringstream oss;
                oss << "Exemplar alignment contains indel not found in the overlap indel set\n"
                    << "\tIndel: " << indelKey
                    << "Exemplar overlap set:\n";
                dump_indel_status(indel_status_map,oss);
                throw blt_exception(oss.str().c_str());
            }
            if (indelKey.isMismatch())
            {
                // this read contains mismatches in the indel buffer
                // we need to recompute candidate alignment
                isRecomputeRequired = true;
            }
            indel_status_map[indelKey].is_present = true;
            indel_status_map[indelKey].isInOriginalAlignment = true;
            validIndels.insert(indelKey);
        }
        if (isRecomputeRequired)
        {
            // Convert the input alignment path to differentiate sequence match (=) and mismatch (X) from
            // 'alignment match' (M). This is required for mismatch processing to work correctly.
            const pos_t readStartPos(unalignedPrefixSize(cal.al.path));
            cal = make_start_pos_alignment(cal.al.pos, readStartPos, cal.al.is_fwd_strand, rseg.read_size(), validIndels);
        }
    }

    // to prevent incompatible alignments, we must put all indels present in the exemplar first in the order list:
    //
    {
        indel_order.clear();
        for (const auto& val : indel_status_map)
        {
            if (val.second.is_present) indel_order.push_back(val.first);
        }
        for (const auto& val : indel_status_map)
        {
            if (! val.second.is_present) indel_order.push_back(val.first);
        }
    }

    // to prevent truncated search, we must put all non-present remove-only indels last in the order list:
    sort_remove_only_indels_last(indel_status_map,indel_order);

#ifdef DEBUG_ALIGN
    std::cerr << "VARMIT exemplar starting indel_status_map:\n";
    dump_indel_status(indel_status_map,std::cerr);

    {
        std::cerr << "VARMIT exemplar starting indel_order:\n";
        const unsigned foo(indel_order.size());
        for (unsigned j(0); j<foo; ++j)
        {
            std::cerr << "no: " << j << " " << indel_order[j] << "\n";
        }
    }
#endif

    // to handle soft-clip and hard-clip in the genomic alignment, we take the
    // soft and hard clip portion off of the alignment before entering
    // the candidate alignment search routine, and add it back in afterwards.
    //
    // TODO -- something less hacky to handle clipping
    //
    unsigned cal_read_length(read_length);
    unsigned hc_lead(0);
    unsigned hc_trail(0);
    unsigned sc_lead(0);
    unsigned sc_trail(0);
    const bool is_input_alignment_clipped(is_clipped(cal.al.path));
    if (is_input_alignment_clipped)
    {
        // exemplar clip condition should only be true for the
        // genomic alignment, which comes first in the exemplar list:
        //
        assert(cal_set.empty());

        apath_clip_clipper(cal.al.path,
                           hc_lead,
                           hc_trail,
                           sc_lead,
                           sc_trail);

        assert(cal_read_length >= (sc_lead+sc_trail));
        cal_read_length-=(sc_lead+sc_trail);
    }

    // launch recursive re-alignment routine starting from the current exemplar alignment:
    static const unsigned startDepth(0);
    static const unsigned startIndelToggleDepth(0);
    static const unsigned startTotalToggleDepth(0);

    candidate_alignment_search(opt, dopt, rseg.getReadIndex(), cal_read_length, indelBuffer,
                               sampleId, realign_buffer_range, cal_set,
                               warn, indel_status_map, haplotypeStatusMap,
                               indel_order, startDepth, startIndelToggleDepth, startTotalToggleDepth,
                               exemplar_pr, opt.max_read_indel_toggle, cal);

    if (is_input_alignment_clipped)
    {
        // un soft-clip candidate alignments:
        std::set<CandidateAlignment> cal_set2(cal_set);
        cal_set.clear();
        for (CandidateAlignment ical : cal_set2)
        {
            apath_clip_adder(ical.al.path,
                             hc_lead,
                             hc_trail,
                             sc_lead,
                             sc_trail);
            cal_set.insert(ical);
        }
    }

    {
        // clear out-of-range alignment candidates:
        std::set<CandidateAlignment> cal_set2(cal_set);
        cal_set.clear();
        for (const CandidateAlignment& ical : cal_set2)
        {
            // check that the alignment is within realign bounds
            if (is_alignment_spanned_by_range(realign_buffer_range,ical.al))
            {
                cal_set.insert(ical);
            }
        }
    }
}



/// standardize input alignment to "unroll" edge insertion and remove edge indels
/// (with exceptions for RNA segment edges adjacent to gap segments)
///
static
alignment
normalizeInputAlignmentIndels(
    const read_segment& rseg)
{
    const std::pair<bool,bool> end_pin(rseg.get_segment_edge_pin());
    const bool is_remove_leading_edge_indels(! end_pin.first);
    const bool is_remove_trailing_edge_indels(! end_pin.second);

    const alignment& inputAlignment(rseg.getInputAlignment());

    const alignment* alignmentPtr(&inputAlignment);
    alignment noEdgeIndelAlignment;
    if ((is_remove_leading_edge_indels || is_remove_trailing_edge_indels) && is_edge_readref_len_segment(inputAlignment.path))
    {
        noEdgeIndelAlignment=matchify_edge_indels(*alignmentPtr,is_remove_leading_edge_indels,is_remove_trailing_edge_indels);
        alignmentPtr=&noEdgeIndelAlignment;
    }

    return *alignmentPtr;
}



void
realignAndScoreRead(
    const starling_base_options& opt,
    const starling_base_deriv_options& dopt,
    const starling_sample_options& sample_opt,
    const reference_contig_segment& ref,
    const known_pos_range& realign_buffer_range,
    const unsigned sampleId,
    read_segment& rseg,
    IndelBuffer& indelBuffer)
{
    if (! rseg.is_valid())
    {
        log_os << "ERROR: invalid alignment path associated with read segment:\n" << rseg;
        exit(EXIT_FAILURE);
    }

    const alignment& inputAlignment(rseg.getInputAlignment());
    assert (! inputAlignment.empty());

    if (! inputAlignment.is_realignable(opt.maxIndelSize)) return;

    if (! check_for_candidate_indel_overlap(realign_buffer_range, rseg, indelBuffer)) return;

    alignment normalizedInputAlignment(normalizeInputAlignmentIndels(rseg));

    const bool isSoftClippedInputAlignment(is_soft_clipped(normalizedInputAlignment.path));
    alignment softClippedInputAlignment;
    if (isSoftClippedInputAlignment)
    {
        softClippedInputAlignment = normalizedInputAlignment;
        normalizedInputAlignment = matchify_edge_soft_clip(softClippedInputAlignment);
    }

    // skip if normalizedInputAlignment has negative start position
    // note this won't come from the read mapper in most cases, but rarely
    // the normalization process could do this
    if (normalizedInputAlignment.pos<0) return;

    // run recursive alignment search starting from normalizedInputAlignment,
    // produce a list of candidate alignments from this search
    //
    // scheme:
    //
    // 1) starting indel set are those overlapped by the input alignment
    //
    // 2) order is defined over the indel set -- indels already present in the
    // input alignment must come first to prevent interference under this
    // scheme
    //
    // 3) when a new candidate alignment overlaps a new indel (ie. not
    // in the starting indel set), that indel is added to the indel set and is
    // pushed onto the end of the indel order.
    //
    // 4) indel status is recorded in the indel_status_map -- toggled
    // on or off to indicate the current state
    //
    // 5) alignments are recorded when recursion depth == indel set
    // size -- note that for this reason the indel off state has to
    // interpreted as "indel not present" rather than "reference", so
    // that all indels can be visited even if some conflict.
    //
    std::set<CandidateAlignment> cal_set;
    mca_warnings warn;

    getCandidateAlignments(opt, dopt, ref, rseg, indelBuffer, sampleId, normalizedInputAlignment,
                           realign_buffer_range, warn, cal_set);

    if ( cal_set.empty() )
    {
        std::ostringstream oss;
        oss << "Empty candidate alignment set while realigning normed input alignment: " << normalizedInputAlignment << "\n";
        throw blt_exception(oss.str().c_str());
    }

    const bool is_incomplete_search(warn.origin_skip || warn.max_toggle_depth);

    // the max_toggle event is too common in genomic resequencing to have a
    // default warning:
    //
    const bool is_max_toggle_warn_enabled(opt.verbosity >= LOG_LEVEL::ALLWARN);
    const bool is_max_toggle_warn(warn.max_toggle_depth && is_max_toggle_warn_enabled);

    if (warn.origin_skip || (is_max_toggle_warn))
    {
        auto writeSkipWarning = [&rseg](
                                    const char* reason)
        {
            log_os << "WARNING: re-alignment skipped some alternate alignments for read: "  << rseg.key()
                   //               << " at buffer_pos: " << (sread.buffer_pos+1) << "\n"
                   << "\treason: " << reason << "\n";
        };

        if (warn.origin_skip) writeSkipWarning("alignments crossed chromosome origin");
        if (is_max_toggle_warn) writeSkipWarning("exceeded max number of indel switches");
    }

    const bool isTestSoftClippedInputAligned(opt.isRetainOptimalSoftClipping && isSoftClippedInputAlignment);
    scoreCandidateAlignmentsAndIndels(opt, dopt, sample_opt, ref,
                                      rseg, indelBuffer, sampleId, cal_set, is_incomplete_search,
                                      isTestSoftClippedInputAligned, softClippedInputAlignment);
}

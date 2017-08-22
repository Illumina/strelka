//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2017 Illumina, Inc.
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

#include "starling_read_util.hh"

#include "blt_util/log.hh"
#include "common/Exceptions.hh"
#include "htsapi/align_path_bam_util.hh"
#include "starling_common/starling_read.hh"

#include <iostream>



void
starling_segmented_read::
pushNewReadSegment(
    const uint16_t size,
    const uint16_t offset,
    const starling_read& sread,
    const alignment& inputAlignment)
{
    _seg_info.emplace_back(size, offset, sread, inputAlignment);
}



read_segment&
starling_segmented_read::
get_segment(const seg_id_t seg_no)
{
    assert(seg_no != 0);
    return _seg_info[seg_no-1];
}



const read_segment&
starling_segmented_read::
get_segment(const seg_id_t seg_no) const
{
    assert(seg_no != 0);
    return _seg_info[seg_no-1];
}



starling_read::
starling_read(
    const bam_record& br,
    const alignment& inputAlignment,
    const MAPLEVEL::index_t inputAlignmentMapLevel)
    : _inputAlignmentMapLevel(inputAlignmentMapLevel),
      _id(0),
      _read_rec(br),
      _full_read(_read_rec.read_size(), 0, *this, inputAlignment)
{
    //
    // additional logic below only applies for spliced RNA-Seq alignments:
    //
    const seg_id_t readSegmentCount(apath_exon_count(inputAlignment.path));
    if (readSegmentCount<=1) return;

    // deal with segmented reads now:
    assert(! is_segmented());
    _segment_ptr.reset(new starling_segmented_read());

    using namespace ALIGNPATH;

    seg_id_t seg_id(1);
    pos_t read_pos(0);
    pos_t ref_pos(inputAlignment.pos);
    pos_t seg_start_read_pos(read_pos);
    pos_t seg_start_ref_pos(ref_pos);
    path_t seg_path;

    const unsigned as(inputAlignment.path.size());
    for (unsigned i(0); i<as; ++i)
    {
        const path_segment& ps(inputAlignment.path[i]);
        const pos_t last_read_pos(read_pos);
        if (is_segment_type_ref_length(ps.type)) ref_pos += ps.length;
        if (is_segment_type_read_length(ps.type)) read_pos += ps.length;

        if (ps.type!=SKIP) seg_path.push_back(ps);

        if (ps.type==SKIP || ((i+1)==as))
        {
            const pos_t end_read_pos( (ps.type==SKIP) ?
                                      last_read_pos : read_pos );
            assert(end_read_pos>seg_start_read_pos);
            const unsigned size(end_read_pos-seg_start_read_pos);
            alignment exonAlignment;
            exonAlignment.path=seg_path;
            exonAlignment.pos=seg_start_ref_pos;
            exonAlignment.is_fwd_strand=inputAlignment.is_fwd_strand;

            _segment_ptr->pushNewReadSegment(size, seg_start_read_pos, *this, exonAlignment);

            seg_id++;
            seg_start_read_pos=read_pos;
            seg_start_ref_pos=ref_pos;
            seg_path.clear();
        }
    }
}



bool
starling_read::
is_tier1_mapping() const
{
    return (MAPLEVEL::TIER1_MAPPED == _inputAlignmentMapLevel);
}



bool
starling_read::
is_tier1or2_mapping() const
{
    return (is_tier1_mapping() ||
            (MAPLEVEL::TIER2_MAPPED == _inputAlignmentMapLevel));
}



uint8_t
starling_read::
map_qual() const
{
    if (!get_full_segment().getInputAlignment().empty())
    {
        return _read_rec.map_qual();
    }
    else
    {
        return 255;
    }
}



void
starling_read::
update_full_segment()
{
    read_segment& fullseg(_full_read);
    assert(! fullseg.is_realigned);

    bool is_new_realignment(false);

    // first determine if anything even needs to be done:
    const unsigned n_seg(segment_count());
    for (unsigned i(0); i<n_seg; ++i)
    {
        const seg_id_t seg_id(i+1);
        const read_segment& rseg(get_segment(seg_id));
        if (! rseg.is_realigned) continue;
        if (rseg.realignment== rseg.getInputAlignment()) continue;
        is_new_realignment=true;
        break;
    }

    if (! is_new_realignment) return;

    // sew together segment reads:
    fullseg.is_realigned=true;
    alignment& fal(fullseg.realignment);
    assert(fal.empty());
    pos_t ref_pos(0);
    for (unsigned i(0); i<n_seg; ++i)
    {
        using namespace ALIGNPATH;
        const seg_id_t seg_id(i+1);
        const read_segment& rseg(get_segment(seg_id));
        const alignment& ral(rseg.is_realigned ? rseg.realignment : rseg.getInputAlignment() );
        assert(! ral.empty());
        if (i==0)
        {
            fal.pos=ral.pos;
            fal.is_fwd_strand=is_fwd_strand();
        }
        else
        {
            // add the exon:
            path_segment ps;
            ps.type = SKIP;
            assert(ral.pos>ref_pos);
            ps.length = ral.pos-ref_pos;
            fal.path.push_back(ps);
        }
        ref_pos=ral.pos;
        for (const path_segment& ps : ral.path)
        {
            if (is_segment_type_ref_length(ps.type)) ref_pos += ps.length;
            fal.path.push_back(ps);
        }
    }
}



void
starling_read::
write_bam(bam_dumper& bamd)
{
    if (is_segmented()) update_full_segment();

    const read_segment& rseg(get_full_segment());

    const alignment& bestAlignment(rseg.getBestAlignment());

    // If the best alignment to exactly the original input alignment, then don't include it in the realigned BAM
    if (bestAlignment == rseg.getInputAlignment()) return;

    // \TODO there should be a soft-clip for the negative position
    // realignment case, for now we skip any bam output for reads where this occurs
    //
    if (bestAlignment.pos < 0) return;

    //
    // write out realigned record:
    //
    bam1_t& br(*_read_rec.get_data());
    bam1_core_t& ca(br.core);

    const bool is_orig_unmapped(_read_rec.is_unmapped());

    // mark mapped bit if necessary:
    if (is_orig_unmapped)
    {
        static const uint8_t unknown_mapq(255);
        ca.flag &= ~(BAM_FLAG::UNMAPPED);
        ca.qual=unknown_mapq;
    }

    // deal with optional fields:
    //

    // if MD or NM records occur (from BWA/novoalign/samtools calmd, etc) it must be
    // taken out before we change the alignment:
    static const char mdtag[] = {'M','D'};
    nuke_bam_aux_field(br,mdtag);
    static const char nmtag[] = {'N','M'};
    nuke_bam_aux_field(br,nmtag);

    // update pos field if it has changed:
    //
    const int32_t orig_pos( is_orig_unmapped ? -1 : ca.pos );
    if (orig_pos != bestAlignment.pos)
    {
        // write current pos to "OP" field if "OP" field does not
        // exist already:
        static const char optag[] = {'O','P'};
        if (nullptr==bam_aux_get(&br,optag))
        {
            assert(orig_pos>=-1);
            const uint32_t out_pos(orig_pos+1);
            bam_aux_append_unsigned(br,optag,out_pos);
        }
        ca.pos=bestAlignment.pos;
    }

    // store orig cigar string if not in aux already (just assume
    // cigar has changed in realignment):
    //
    static const char octag[] = {'O','C'};
    if ((! is_orig_unmapped) && (nullptr==bam_aux_get(&br,octag)))
    {
        std::string _oc_cigar;
        apath_to_cigar(rseg.getInputAlignment().path,_oc_cigar);
        bam_aux_append(&br,octag,'Z', (_oc_cigar.size()+1),(uint8_t*) (_oc_cigar.c_str()));
    }

    // update cigar field:
    edit_bam_cigar(bestAlignment.path,br);

    bam_update_bin(br);
    bamd.put_record(&br);
}



std::ostream&
operator<<(std::ostream& os,
           const starling_read& sr)
{
    os << "STARLING_READ id: " << sr.id()
       << " mapping_type: " << MAPLEVEL::get_label(sr.getInputAlignmentMapLevel())
       << "\n";

    os << "full_segment_info:\n";
    os << sr.get_full_segment();

    const seg_id_t sc(sr.segment_count());
    for (unsigned i(0); i<sc; ++i)
    {
        const seg_id_t seg_id(i+1);
        os << "partial_segment " << static_cast<unsigned>(seg_id) << " :\n";
        short_report(os,sr.get_segment(seg_id));
    }

    return os;
}

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

#include "starling_read_util.hh"

#include "blt_util/log.hh"
#include "common/Exceptions.hh"
#include "htsapi/align_path_bam_util.hh"
#include "starling_common/starling_read.hh"

#include <iostream>



starling_read::
starling_read(
    const bam_record& br,
    const alignment& inputAlignment,
    const MAPLEVEL::index_t inputAlignmentMapLevel,
    const align_id_t readIndex)
    : _inputAlignmentMapLevel(inputAlignmentMapLevel),
      _readIndex(readIndex),
      _read_rec(br),
      _full_read(_read_rec.read_size(), 0, *this, inputAlignment)
{
    const seg_id_t exonCount(apath_exon_count(inputAlignment.path));
    if (exonCount <= 1) return;

    //
    // additional logic below only applies for spliced RNA-Seq alignments:
    //

    using namespace ALIGNPATH;

    pos_t readPos(0);
    pos_t refPos(inputAlignment.pos);

    pos_t exonStartReadPos(readPos);
    pos_t exonStartRefPos(refPos);
    path_t exonAlignmentPath;

    const unsigned as(inputAlignment.path.size());
    for (unsigned i(0); i<as; ++i)
    {
        const path_segment& ps(inputAlignment.path[i]);
        const pos_t lastReadPos(readPos);
        if (is_segment_type_ref_length(ps.type)) refPos += ps.length;
        if (is_segment_type_read_length(ps.type)) readPos += ps.length;

        if (ps.type!=SKIP) exonAlignmentPath.push_back(ps);

        if (ps.type==SKIP || ((i+1)==as))
        {
            const pos_t endReadPos( (ps.type==SKIP) ?
                                    lastReadPos : readPos );
            assert(endReadPos>exonStartReadPos);
            const unsigned exonReadSize(endReadPos-exonStartReadPos);
            alignment exonAlignment;
            exonAlignment.path=exonAlignmentPath;
            exonAlignment.pos=exonStartRefPos;
            exonAlignment.is_fwd_strand=inputAlignment.is_fwd_strand;

            _exonInfo.emplace_back(exonReadSize, exonStartReadPos, *this, exonAlignment);

            exonStartReadPos=readPos;
            exonStartRefPos=refPos;
            exonAlignmentPath.clear();
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



void
starling_read::
update_full_segment()
{
    read_segment& fullseg(_full_read);
    assert(! fullseg.is_realigned);

    bool is_new_realignment(false);

    // first determine if anything even needs to be done:
    const unsigned exonCount(getExonCount());
    for (unsigned exonIndex(0); exonIndex<exonCount; ++exonIndex)
    {
        const seg_id_t segmentIndex(exonIndex+1);
        const read_segment& rseg(get_segment(segmentIndex));
        if (! rseg.is_realigned) continue;
        if (rseg.realignment == rseg.getInputAlignment()) continue;
        is_new_realignment=true;
        break;
    }

    if (! is_new_realignment) return;

    // sew together segment reads:
    fullseg.is_realigned=true;
    alignment& fullRealignment(fullseg.realignment);
    assert(fullRealignment.empty());
    pos_t refPos(0);
    for (unsigned exonIndex(0); exonIndex<exonCount; ++exonIndex)
    {
        using namespace ALIGNPATH;
        const seg_id_t segmentIndex(exonIndex+1);
        const read_segment& exonSegment(get_segment(segmentIndex));
        const alignment& exonAlignment(exonSegment.getBestAlignment());
        assert(! exonAlignment.empty());
        if (exonIndex==0)
        {
            fullRealignment.pos=exonAlignment.pos;
            fullRealignment.is_fwd_strand=is_fwd_strand();
        }
        else
        {
            // add the exon:
            path_segment ps;
            ps.type = SKIP;
            assert(exonAlignment.pos>refPos);
            ps.length = exonAlignment.pos-refPos;
            fullRealignment.path.push_back(ps);
        }
        refPos=exonAlignment.pos;
        for (const path_segment& ps : exonAlignment.path)
        {
            if (is_segment_type_ref_length(ps.type)) refPos += ps.length;
            fullRealignment.path.push_back(ps);
        }
    }
}



void
starling_read::
write_bam(bam_dumper& bamd)
{
    if (isSpliced()) update_full_segment();

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

        // Do lots of ugly casting on svStr to fit htsapi signature. Usage is actually const in htslib:
        bam_aux_append(&br,octag,'Z', (_oc_cigar.size()+1),
                       reinterpret_cast<uint8_t*>(const_cast<char*>(_oc_cigar.c_str())));
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
    os << "STARLING_READ id: " << sr.getReadIndex()
       << " mapping_type: " << MAPLEVEL::get_label(sr.getInputAlignmentMapLevel())
       << "\n";

    os << "full_segment_info:\n";
    os << sr.get_full_segment();

    const seg_id_t exonCount(sr.getExonCount());
    for (unsigned exonIndex(0); exonIndex<exonCount; ++exonIndex)
    {
        const seg_id_t segmentIndex(exonIndex+1);
        os << "partial_segment " << static_cast<unsigned>(segmentIndex) << " :\n";
        short_report(os,sr.get_segment(segmentIndex));
    }

    return os;
}

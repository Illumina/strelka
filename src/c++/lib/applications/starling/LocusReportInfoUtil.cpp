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

///
/// \author Chris Saunders
///

#include "LocusReportInfoUtil.hh"
#include "common/Exceptions.hh"

#include <sstream>



/// helper for reference_contig_segment
static
void
append_ref_subseq(
    const reference_contig_segment& ref,
    const pos_t start_pos,
    const pos_t end_pos,
    std::string& out_seq)
{
    for (pos_t p(start_pos); p<end_pos; ++p)
    {
        out_seq += ref.get_base(p);
    }
}



/// helper for reference_contig_segment
static
void
copy_ref_subseq(
    const reference_contig_segment& ref,
    const pos_t start_pos,
    const pos_t end_pos,
    std::string& out_seq)
{
    out_seq.clear();
    append_ref_subseq(ref, start_pos, end_pos, out_seq);
}



void
getLocusReportInfoFromAlleles(
    const reference_contig_segment& ref,
    const std::vector<GermlineIndelAlleleInfo>& indelAlleles,
    unsigned commonPrefixLength,
    OrthogonalAlleleSetLocusReportInfo& locusReportInfo)
{
    assert(indelAlleles.size()>0);

    locusReportInfo.altAlleles.resize(indelAlleles.size());
    std::string& vcfRefSeq(locusReportInfo.vcfRefSeq);

    // overlapping breakpoints are not allowed, but we do handle singles:
    if ((indelAlleles.size()==1) and (indelAlleles[0].indelKey.is_breakpoint()))
    {
        const IndelKey& indelKey(indelAlleles[0].indelKey);
        std::string& vcfAltSeq(locusReportInfo.altAlleles[0].vcfAltSeq);

        if       (indelKey.type == INDEL::BP_LEFT)
        {
            copy_ref_subseq(ref, indelKey.pos-1, indelKey.pos, vcfRefSeq);
            vcfAltSeq = vcfRefSeq + indelAlleles[0].breakpointInsertSeq + '.';
            locusReportInfo.vcfPos = indelKey.pos;
        }
        else if (indelKey.type == INDEL::BP_RIGHT)
        {
            copy_ref_subseq(ref, indelKey.pos, indelKey.pos+1, vcfRefSeq);
            vcfAltSeq = '.' + indelAlleles[0].breakpointInsertSeq + vcfRefSeq;
            locusReportInfo.vcfPos = indelKey.pos+1;
        }
        else
        {
            assert(false and "Breakpoint allele type expected");
        }
    }
    else
    {
        bool isFirst(true);
        known_pos_range2 alleleSetRange;
        for (const auto& indelAllele : indelAlleles)
        {
            const IndelKey& indelKey(indelAllele.indelKey);
            assert(not indelKey.is_breakpoint());

            if (isFirst)
            {
                alleleSetRange.set_range(indelKey.pos, indelKey.right_pos());
                isFirst = false;
            }
            else
            {
                alleleSetRange.merge_range(known_pos_range2(indelKey.pos, indelKey.right_pos()));
            }
        }

        alleleSetRange.set_begin_pos(alleleSetRange.begin_pos());

        // minus 1 for vcf preceding base, and plus 1 for 0-index to 1-index shift equals -- same value:
        // shift the pos if commonPrefixLength > 1
        locusReportInfo.vcfPos = alleleSetRange.begin_pos() + commonPrefixLength;

        copy_ref_subseq(ref, alleleSetRange.begin_pos() - 1 + commonPrefixLength, alleleSetRange.end_pos(), vcfRefSeq);

        const unsigned altAlleleCount(indelAlleles.size());
        for (unsigned altAlleleIndex(0); altAlleleIndex < altAlleleCount; ++altAlleleIndex)
        {
            const IndelKey& indelKey(indelAlleles[altAlleleIndex].indelKey);
            std::string& vcfAltSeq(locusReportInfo.altAlleles[altAlleleIndex].vcfAltSeq);

            copy_ref_subseq(ref, alleleSetRange.begin_pos() - 1, indelKey.pos, vcfAltSeq);
            unsigned leadingPad(indelKey.pos - (alleleSetRange.begin_pos() - 1));

            vcfAltSeq += indelKey.insert_seq();

            append_ref_subseq(ref, indelKey.right_pos(), alleleSetRange.end_pos(), vcfAltSeq);
            unsigned trailingPad(alleleSetRange.end_pos() - indelKey.right_pos());

            // trim the common prefix
            vcfAltSeq = vcfAltSeq.substr(commonPrefixLength);

            // get corresponding CIGAR string:
            auto& vcfCigar(locusReportInfo.altAlleles[altAlleleIndex].vcfCigar);
            setIndelAlleleCigar(leadingPad, trailingPad, commonPrefixLength, indelKey, vcfCigar);
        }

        // compare all ALT sequences to sanity check that they don't match
        /// TODO how long do we need to keep this assertion here?
        std::set<std::string> vcfAltSeqCheck;
        for (unsigned altAlleleIndex(0); altAlleleIndex < altAlleleCount; ++altAlleleIndex)
        {
            const std::string& vcfAltSeq(locusReportInfo.altAlleles[altAlleleIndex].vcfAltSeq);
            const auto insertRetval = vcfAltSeqCheck.insert(vcfAltSeq);
            if (not insertRetval.second)
            {
                using namespace illumina::common;

                std::ostringstream oss;
                oss << "Repeated VCF ALT value in indel record: '" << vcfAltSeq << "'\n";
                oss << "\tIndel Alleles (" << indelAlleles.size() << "):\n";
                for (const auto& indelAllele : indelAlleles)
                {
                    oss << indelAllele << "\n";
                }
                oss << "\n";
                BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
            }
        }
    }
}

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

#include "vcf_record_util.hh"
#include "common/Exceptions.hh"

#include <sstream>


bool
isExpectedVcfReference(
    const reference_contig_segment& ref,
    const vcf_record& vcfRecord)
{
    const unsigned vcfRefSize(vcfRecord.ref.size());
    for (unsigned vcfRefIndex(0); vcfRefIndex < vcfRefSize; ++vcfRefIndex)
    {
        const char fastaRefBase(ref.get_base(vcfRecord.pos - 1 + vcfRefIndex));
        const char vcfRefBase(vcfRecord.ref[vcfRefIndex]);
        if (fastaRefBase == vcfRefBase) continue;
        if ((fastaRefBase == 'N') || (vcfRefBase == 'N')) continue;
        return false;
    }

    return true;
}



void
assertExpectedVcfReference(
    const reference_contig_segment& ref,
    const vcf_streamer& vcfStreamer)
{
    const vcf_record* vptr(vcfStreamer.get_record_ptr());
    assert(vptr != nullptr);
    const vcf_record& vcfRecord(*vptr);

    if (isExpectedVcfReference(ref,vcfRecord)) return;

    std::string fastaReferenceSegment;
    ref.get_substring(vcfRecord.pos-1, vcfRecord.ref.size(), fastaReferenceSegment);

    std::ostringstream oss;
    oss << "Input VCF record REF value is not compatible with genome reference:\n";
    vcfStreamer.report_state(oss);
    oss << "Genome reference: '" << fastaReferenceSegment << "'\n";
    oss << "VCF record REF value:    '" << vcfRecord.ref << "'\n";
    oss << "Please ensure that the input VCF comes from the appropriate reference genome\n";
    BOOST_THROW_EXCEPTION(illumina::common::GeneralException(oss.str()));
}

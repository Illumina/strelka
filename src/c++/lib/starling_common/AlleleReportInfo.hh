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

#pragma once

#include "blt_common/MapqTracker.hh"
#include "blt_util/fastRanksum.hh"
#include "blt_util/MeanTracker.hh"
#include "starling_common/IndelKey.hh"

#include <iosfwd>
#include <string>


namespace SimplifiedIndelReportType
{

/// some components of the indel reporting need to reduce all alternate alleles to the following
/// simplified states:
enum index_t
{
    INSERT,
    DELETE,
    SWAP,
    BREAKPOINT,
    OTHER
};

inline
index_t
getRateType(
    const IndelKey& indelKey)
{
    if     (indelKey.isPrimitiveDeletionAllele())
    {
        return DELETE;
    }
    else if (indelKey.isPrimitiveInsertionAllele())
    {
        return INSERT;
    }
    else if (indelKey.type == INDEL::INDEL)
    {
        return SWAP;
    }
    else if (indelKey.is_breakpoint())
    {
        return BREAKPOINT;
    }
    else
    {
        return OTHER;
    }
}
}


/// \brief Allele summary information which is shared between all samples.
///
struct AlleleReportInfo
{
    AlleleReportInfo() {}

    bool
    isRepeatUnit() const
    {
        return (! repeatUnit.empty());
    }

    void dump(std::ostream& os) const;

    std::string repeatUnit;
    unsigned repeatUnitLength = 0;
    unsigned refRepeatCount = 0;
    unsigned indelRepeatCount = 0;
    unsigned interruptedHomopolymerLength = 0;
    unsigned contextCompressability = 0; /// max adjacent sequence length encodable by 5 Ziv-Lempel-77 keywords

    // not directly reported, but handy to have pre-calculated:
    SimplifiedIndelReportType::index_t it = SimplifiedIndelReportType::OTHER;
};

std::ostream& operator<<(std::ostream& os, const AlleleReportInfo& obj);


/// \brief Allele summary information which is specific to each sample.
///
struct AlleleSampleReportInfo
{
    /// TODO STREL-125 sample_report_info is still designed for only one alt allele

    AlleleSampleReportInfo() {}

    unsigned n_confident_ref_reads = 0;
    unsigned n_confident_indel_reads = 0;
    unsigned n_confident_alt_reads = 0;

    /// Number of reads which cannot be confidently assigned to one of the modeled alleles or are otherwise low quality
    unsigned n_other_reads = 0;

    /// Indel depth estimated from the pileup depth at the position preceding the indel
    unsigned indelLocusDepth = 0;

    // same as above, but by strand
    unsigned n_confident_ref_reads_fwd = 0;
    unsigned n_confident_indel_reads_fwd = 0;
    unsigned n_confident_alt_reads_fwd = 0;
    unsigned n_other_reads_fwd = 0;
    unsigned n_confident_ref_reads_rev = 0;
    unsigned n_confident_indel_reads_rev = 0;
    unsigned n_confident_alt_reads_rev = 0;
    unsigned n_other_reads_rev = 0;

    MapqTracker mapqTracker;

    fastRanksum readpos_ranksum;

    /// supports RNA-Seq EVS feature
    MeanTracker distanceFromReadEdge;

    unsigned total_confident_reads() const
    {
        return n_confident_alt_reads + n_confident_indel_reads + n_confident_ref_reads;
    }

    /// Debugging printer
    void dump(std::ostream& os) const;
};

/// Debugging printer
std::ostream& operator<<(std::ostream& os, const AlleleSampleReportInfo& obj);

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
/// \author Mitch Bekritsky
///

#pragma once

#include "blt_util/known_pos_range2.hh"
#include "blt_util/blt_types.hh"

#include "htsapi/vcf_record.hh"

#include "starling_common/starling_diploid_indel.hh"
#include "starling_common/IndelKey.hh"

#include "boost/icl/discrete_interval.hpp"
#include "boost/icl/interval_map.hpp"


#include <iosfwd>
#include <set>


namespace GENOTYPE_STATUS
{
enum genotype_t
{
    UNKNOWN,
    HOMREF,
    HET,
    HETALT,
    HOMALT,
    SIZE
};

inline
const char*
label(
    const genotype_t gt)
{
    switch (gt)
    {
    case UNKNOWN:
        return "Unknown";
    case HOMREF:
        return "Homref";
    case HET:
        return "Het";
    case HETALT:
        return "Hetalt";
    case HOMALT:
        return "Homalt";
    default:
        return "xxx";
    }
}
}


struct IndelGenotype
{
    // alts are stored in a vector to maintain the ordering
    // of alts in the VCF file
    typedef std::vector<IndelKey> alt_vec_t;
    IndelGenotype(
        const vcf_record& in_vcfr)
        : vcfr(in_vcfr)
    {
        pos = in_vcfr.pos;
        extractAlts();
        assignGenotype();
        assignFilter();
    }

    bool
    operator==(const IndelGenotype& rhs) const
    {
        if (vcfr.pos == rhs.vcfr.pos && genotype == rhs.genotype) return true;
        return false;
    }

    // return true if the vcf record contains a match to indelKey
    bool
    altMatch(
        const IndelKey& indelKey) const;

    vcf_record vcfr;
    pos_t pos;
    alt_vec_t alts;
    unsigned max_delete_length;
    std::string gt_string;
    std::string filter;
    GENOTYPE_STATUS::genotype_t genotype = GENOTYPE_STATUS::HOMREF;

private:
    void extractAlts();
    void assignGenotype();
    void assignFilter();
};


struct IndelGenotypeSort
{
    bool
    operator()(
        const IndelGenotype& lhs,
        const IndelGenotype& rhs)
    {

        if (lhs.vcfr.pos < rhs.vcfr.pos) return true;
        if (lhs.vcfr.pos == rhs.vcfr.pos)
        {
            // respect GENOTYPE_STATUS genotype ordering
            if (lhs.genotype < rhs.genotype) return true;
        }

        // I'm assuming that VCF records beginning at the same position
        // do not have the same genotype, i.e. there's no such thing as
        // two records starting at the same site with the same genotype
        return false;
    }
};


struct RecordTracker
{
    typedef std::set<IndelGenotype, IndelGenotypeSort> indel_value_t;

    bool
    empty() const
    {
        return _records.empty();
    }

    void
    clear()
    {
        _records.clear();
    }

    bool
    intersectingRecord(
        const pos_t pos,
        indel_value_t& records) const
    {
        return intersectingRecordImpl(pos, pos + 1, records);
    }

    bool
    intersectingRecord(
        const known_pos_range2 range,
        indel_value_t& records) const
    {
        return intersectingRecordImpl(range.begin_pos(), range.end_pos(), records);
    }

    void
    addVcfRecord(
        const vcf_record& vcfRecord);

    void
    dump(
        std::ostream& os) const;

    unsigned
    size() const
    {
        return _records.size();
    }

private:
    typedef boost::icl::interval_map<pos_t, indel_value_t> indel_record_t;
    typedef boost::icl::discrete_interval<pos_t> interval_t;

    bool
    intersectingRecordImpl(
        const pos_t beginPos,
        const pos_t endPos,
        RecordTracker::indel_value_t& records) const;

    indel_record_t _records;
};

std::ostream&
operator<<(
    std::ostream& os,
    const IndelGenotype& gt);

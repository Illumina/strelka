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
/// \author Mitch Bekritsky
///

#include "blt_util/RecordTracker.hh"
#include <iostream>

std::ostream&
operator<<(
    std::ostream& os,
    const IndelVariant& var)
{
    // os << "Original ref string: " << var.ref_string << "\n";
    // os << "Original alt string: " << var.alt_string << "\n";
    os << "Indel variant type: " << INDEL::get_index_label(var.type) << "\n";
    os << "Indel length: " << var.length << "\n";
    if (var.type == INDEL::INSERT)
    {
        os << "Inserted sequence: " << var.insert_sequence << "\n";
    }

    return os;
}

std::ostream&
operator<<(
    std::ostream& os,
    const IndelGenotype& gt)
{
    os << "VCF record: " << gt.vcfr;
    os << "Genotype: " << STAR_DIINDEL::label(gt.genotype) << "\n";
    os << "Observed alts:\n";

    unsigned alt_index(1);
    for (const auto& alt : gt.alts)
    {
        os << "Alt " << alt_index << ":\n";
        os << alt;
        alt_index++;
    }

    return os;
}

bool
IndelGenotype::
altMatch(
    const indel_key& ik,
    const indel_data& id) const
{
    for (const auto alt : alts)
    {
        if (ik.type == alt.type && ik.length == alt.length)
        {
            if (ik.type == INDEL::DELETE) return true;
            if (id.get_insert_seq() == alt.insert_sequence) return true;
        }
    }
    return false;
}


bool
RecordTracker::
intersectingRecordImpl(
    const pos_t beginPos,
    const pos_t endPos,
    RecordTracker::indel_value_t& records) const
{
    if (_records.empty()) return false;

    interval_t search_interval = boost::icl::construct<interval_t>(beginPos, endPos, boost::icl::interval_bounds::right_open());

    auto resultIter  = _records.find(search_interval);
    if (resultIter == _records.end()) return false;

    // std::cout << "Search interval overlaps known variant\n";
    // std::cout << "Search interval: " << search_interval << std::endl;
    // std::cout << "Known variant interval: " << resultIter->first << std::endl;
    // std::cout << "Variants:\n";
    // for(const auto& var : resultIter->second)
    // {
    // 	std::cout << "\t" << var << std::endl;
    // }

    records = resultIter->second;

    return true;
}

void
RecordTracker::
addVcfRecord(
    const vcf_record& vcfRecord)
{
    // overlapping vcf records are allowed

    pos_t start_pos = vcfRecord.pos;
    pos_t end_pos   = start_pos + 1;
    // check to see if REF or ALT field has length > 1 (i.e. record is indel)

    unsigned ref_length = vcfRecord.ref.size();

    // currently only used to add VCFs with single alts
    assert(vcfRecord.alt.size() == 1 && vcfRecord.is_indel());

    // if it's an indel, populate the IndelVariant struct
    unsigned alt_length = vcfRecord.alt.begin()->size();

    assert(alt_length != ref_length);
    // AFAIK, the only way for alt_length and ref_length to be
    // equal would be for the record to be a SNP

    IndelVariant variant_instance;
    IndelGenotype genotype_instance;
    genotype_instance.vcfr = vcfRecord;

    variant_instance.ref_string = vcfRecord.ref;
    variant_instance.alt_string = *(vcfRecord.alt.begin());

    if (ref_length < alt_length)
    {
        variant_instance.type = INDEL::INSERT;
        variant_instance.length = alt_length - ref_length;
        variant_instance.insert_sequence = vcfRecord.alt.begin()->substr(1);
    }
    else
    {
        variant_instance.type = INDEL::DELETE;
        variant_instance.length = ref_length - alt_length;
    }

    genotype_instance.alts.insert(variant_instance);

    // deletion has length >1 in record tracker
    if (ref_length > 0 && alt_length == 1)
    {
        end_pos = start_pos + ref_length - 1;
    }
    // insertions and SNPs do not have length >1

    interval_t interval = boost::icl::discrete_interval<pos_t>::right_open(start_pos, end_pos);
    indel_value_t record_set;
    record_set.insert(genotype_instance);

    _records.insert(std::make_pair(interval, record_set));
}

void
RecordTracker::
dump(
    std::ostream& os) const
{
    os << "RecordTracker\n";
    for (const auto& record : _records)
    {
        os << record.first << " => {";
        for (const auto& gt : record.second)
        {
            os << "{" << gt << "};";
        }
        os << "}\n";
    }
}

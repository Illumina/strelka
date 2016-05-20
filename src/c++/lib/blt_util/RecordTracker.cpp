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


// extract an entry from the format field in a VCF record
// and returns the value in the first VCF sample (that's
// easy  to change if needed)
void
extractFormatEntry(
    const std::string& vcf_line,
    const std::string& field,
    std::string& value,
    bool strict = true)
{
    static const char field_sep('\t');
    static const char format_sep(':');
    static const unsigned format_index(8);

    size_t field_start(0);
    size_t next_start(0);

    // fast-forward to format field
    for (unsigned i(0); i < format_index; ++i)
    {
        field_start = vcf_line.find_first_of(field_sep, field_start) + 1;
    }
    next_start = vcf_line.find_first_of(field_sep, field_start) + 1;
    std::string format_string = vcf_line.substr(field_start, next_start - field_start - 1);

    // go to first sample field
    field_start = next_start;
    next_start  = vcf_line.find_first_of(field_sep, field_start);
    std::string sample_string = vcf_line.substr(field_start, next_start - field_start - 1);
    // N.B. if next_start is std::string::npos (i.e. only one sample in VCF), this will grab til
    // the end of the line

    size_t entry_start (format_string.find(field));

    if (strict) assert(entry_start != std::string::npos);

    // count how many separators you need to pass before getting to the right field
    size_t entry_index (std::count(format_string.begin(), format_string.begin() + entry_start, format_sep));

    size_t sample_entry_pos(0);
    for (unsigned i(0); i < entry_index; ++i)
    {
        sample_entry_pos = sample_string.find_first_of(format_sep, sample_entry_pos) + 1;
    }
    size_t next_pos = sample_string.find_first_of(format_sep, sample_entry_pos) + 1;

    value = sample_string.substr(sample_entry_pos, next_pos - sample_entry_pos - 1);
}


// extract genotype from VCF entry
std::string
getGenotypeString(
    const std::string& vcf_line)
{
    std::string genotype;
    extractFormatEntry(vcf_line, "GT", genotype);
    return genotype;
}


IndelVariant::
IndelVariant(
    const vcf_record& vcfr,
    const std::string& alt_instance)
{
    unsigned ref_length = vcfr.ref.size();
    unsigned alt_length = alt_instance.size();

    assert(alt_length != ref_length);
    // AFAIK, the only way for alt_length and ref_length to be
    // equal would be for the record to be a SNP

    alt_string = alt_instance;
    ref_string = vcfr.ref;

    if (ref_length < alt_length)
    {
        type = INDEL::INSERT;
        length = alt_length - ref_length;
        insert_sequence = alt_instance.substr(1);
    }
    else
    {
        type = INDEL::DELETE;
        length = ref_length - alt_length;
    }
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


void
IndelGenotype::
extractAlts()
{
    max_delete_length = 0;
    for(const auto& alt_str : vcfr.alt)
    {
        IndelVariant variant_instance = IndelVariant(vcfr, alt_str);
        alts.push_back(variant_instance);

        if(variant_instance.type == INDEL::DELETE)
        {
            max_delete_length = std::max(max_delete_length, variant_instance.length);
        }
    }
}


void
IndelGenotype::
assignGenotype()
{
    gt_string = getGenotypeString(vcfr.line);
    assert(!gt_string.empty());

    size_t gt_string_delim_pos = gt_string.find_first_of("|/");

    assert(gt_string_delim_pos != std::string::npos);

    int allele_1 = stoi(gt_string.substr(0, gt_string_delim_pos));
    int allele_2 = stoi(gt_string.substr(gt_string_delim_pos + 1));

    // no negative genotype indices
    assert (allele_1 >= 0 && allele_2 >= 0);

    // alt indices must be less than or equal to the the total number of alts
    assert (allele_1 <= (int) alts.size() && allele_2 <= (int) alts.size());

    unsigned hom_count = (allele_1 > 0) + (allele_2 > 0);

    if (hom_count == 0) // i.e. 0/0
    {
        genotype = GENOTYPE_STATUS::HOMREF;
    }
    else if (hom_count == 1) // i.e. 0/1, 0/2, 1/0, etc.
    {
        genotype = GENOTYPE_STATUS::HET;
    }
    else
    {
        if (allele_1 == allele_2) // i.e. 1/1
        {
            genotype = GENOTYPE_STATUS::HOMALT;
        }
        else // i.e. 1/2
        {
            genotype = GENOTYPE_STATUS::HETALT;
        }
    }
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
    // currently only used to add indels with single alts
    assert(vcfRecord.alt.size() == 1 && vcfRecord.is_indel());

    IndelGenotype genotype_instance(vcfRecord);

    // since these intervals are right-open, the end position of the interval
    // is one past the final reference base in the indel
    pos_t start_pos = vcfRecord.pos;
    pos_t end_pos   = start_pos + genotype_instance.max_delete_length + 1;

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
    os << "Pos: " << gt.pos << "\n";
    os << "Genotype: " << GENOTYPE_STATUS::label(gt.genotype) << "\n";
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

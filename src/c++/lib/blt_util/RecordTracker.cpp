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
/// \author Mitch Bekritsky
///

#include "blt_util/RecordTracker.hh"
#include <iostream>
#include <array>



static
size_t
extractVcfField(
    const std::string& vcf_line,
    unsigned field_index,
    std::string& field,
    size_t start_from = 0)
{
    static const char field_sep('\t');

    size_t field_start = start_from;

    // fast-forward to field_index
    for (unsigned i(0); i < field_index; ++i)
    {
        field_start = vcf_line.find_first_of(field_sep, field_start + 1);
    }

    if (field_start != std::string::npos)
    {
        size_t next_start(vcf_line.find_first_of(field_sep, field_start + 1));
        field = vcf_line.substr(field_start + 1, next_start - field_start - 1);

        return next_start + 1;
    }

    // if the next occurence of the tab is std::string::npos, then this field is the
    // last field in the record
    field = vcf_line.substr(start_from);

    return std::string::npos;
}

// extract an entry from the format field in a VCF record
// and returns the value in the first VCF sample (that's
// easy  to change if needed)
static
void
extractFormatEntry(
    const std::string& vcf_line,
    const std::string& field,
    std::string& value,
    bool strict = true)
{
    static const char format_sep(':');
    static const unsigned format_index(8);

    std::string format_string, sample_string;
    size_t formatEnd(extractVcfField(vcf_line, format_index, format_string));
    extractVcfField(vcf_line, 1, sample_string, formatEnd);

    // VCF record must have a sample entry
    assert(!sample_string.empty());

    size_t entry_start (format_string.find(field));
    // if strict is true, require that field is a key in the FORMAT string
    if (strict)
    {
        assert(entry_start != std::string::npos);
    }

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
static
std::string
getGenotypeString(
    const std::string& vcf_line)
{
    std::string genotype;
    extractFormatEntry(vcf_line, "GT", genotype);
    return genotype;
}



static
bool
testIsFirstBaseMatching(
    const std::string& a,
    const std::string& b)
{
    if (a.empty() or b.empty()) return false;
    return (a[0] == b[0]);
}



static
IndelKey
convertVcfAltToIndel(
    const vcf_record& vcfr,
    const std::string& alt_instance)
{
    unsigned ref_length = vcfr.ref.size();
    unsigned alt_length = alt_instance.size();

    assert(alt_length != ref_length);
    // AFAIK, the only way for alt_length and ref_length to be
    // equal would be for the record to be a SNP

    const bool isFirstBaseMatching(testIsFirstBaseMatching(vcfr.ref, alt_instance));
    const pos_t startPosOffset((isFirstBaseMatching ? 1 : 0));
    IndelKey indelKey;
    indelKey.pos = (vcfr.pos -1) + startPosOffset;
    indelKey.type = INDEL::INDEL;
    indelKey.deletionLength = ref_length - startPosOffset;
    indelKey.insertSequence = alt_instance.substr(startPosOffset);

    return indelKey;
}



bool
IndelGenotype::
altMatch(
    const IndelKey& indelKey) const
{
    for (const auto& alt : alts)
    {
        if (indelKey == alt) return true;
    }
    return false;
}


void
IndelGenotype::
extractAlts()
{
    max_delete_length = 0;
    for (const auto& alt_str : vcfr.alt)
    {
        alts.push_back(convertVcfAltToIndel(vcfr, alt_str));
        max_delete_length = std::max(max_delete_length, alts.back().deletionLength);
    }
}


void
IndelGenotype::
assignGenotype()
{
    gt_string = getGenotypeString(vcfr.line);
    assert(!gt_string.empty());

    // we could move this into the array loop, but I want to keep
    // the assert for now, which is easier to have out here
    size_t gt_string_delim_pos = gt_string.find_first_of("|/");
    assert(gt_string_delim_pos != std::string::npos);

    std::array <int, 2> alleles;
    size_t start_pos = 0;
    unsigned hom_count = 0;

    // parse alleles
    for (unsigned i = 0; i < alleles.size(); ++i)
    {
        std::string allele_instance(gt_string.substr(start_pos, gt_string_delim_pos));

        // currently assuming that alleles can either be an integer or a '.',
        // nothing else
        if (allele_instance != ".")
        {
            alleles[i] = stoi(allele_instance);
        }
        else
        {
            alleles[i] = 0;
        }
        start_pos += gt_string_delim_pos + 1;

        // no negative genotype indices and alt indices must be less than or
        // equal to the the total number of alts
        assert(alleles[i] >= 0 && alleles[i] <= (int) alts.size());
        hom_count += (alleles[i] == 0);
    }

    // assign genotype status (homref, het, homalt, hetalt)
    if (hom_count == 2) // i.e. 0/0
    {
        genotype = GENOTYPE_STATUS::HOMREF;
    }
    else if (hom_count == 1) // i.e. 0/1, 0/2, 1/0, etc.
    {
        genotype = GENOTYPE_STATUS::HET;
    }
    else
    {
        if (alleles[0] == alleles[1]) // i.e. 1/1
        {
            genotype = GENOTYPE_STATUS::HOMALT;
        }
        else // i.e. 1/2
        {
            genotype = GENOTYPE_STATUS::HETALT;
        }
    }
}

void
IndelGenotype::
assignFilter()
{
    static const unsigned filter_index = 6;
    extractVcfField(vcfr.line, filter_index, filter);
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

    // exclude non-passing variants here for now (may want to add all
    // records and leave treatment as downstream decision)
    if (genotype_instance.filter == "PASS" || genotype_instance.filter == ".")
    {
        // since these intervals are right-open, the end position of the interval
        // is one past the final reference base in the indel
        pos_t start_pos = vcfRecord.pos;
        pos_t end_pos   = start_pos + genotype_instance.max_delete_length + 1;

        interval_t interval = interval_t::right_open(start_pos, end_pos);
        indel_value_t record_set;
        record_set.insert(genotype_instance);

        _records.insert(std::make_pair(interval, record_set));
    }
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
    const IndelGenotype& gt)
{
    os << "Pos: " << gt.pos << "\n";
    os << "Genotype: " << GENOTYPE_STATUS::label(gt.genotype) << "\n";
    os << "Filter: " << gt.filter << "\n";
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

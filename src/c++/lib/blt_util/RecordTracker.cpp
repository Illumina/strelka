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

bool
RecordTracker::
intersectingRecordImpl(
    const pos_t beginPos,
    const pos_t endPos,
    std::set<std::string>& records) const
{
	if (_records.empty()) return false;

	interval_t search_interval = boost::icl::construct<interval_t>(beginPos, endPos, boost::icl::interval_bounds::right_open());

	auto resultIter  = _records.find(search_interval);
	if (resultIter == _records.end()) return false;

	// std::cout << "Search interval overlaps known variant\n";
	// std::cout << "Search interval: " << search_interval << std::endl;
	// std::cout << "Known variant interval: " << resultIter->first << std::endl;
	// std::cout << "Variants:\n";
	// for(value_t::const_iterator it = resultIter->second.begin(); it != resultIter->second.end(); ++it)
	// {
	// 	std::cout << "\t" << *it << std::endl;
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
	assert(vcfRecord.alt.size() == 1);
	unsigned alt_length = vcfRecord.alt.begin()->size();

	// deletion has length >1 in record tracker
	if (ref_length > 0 && alt_length == 1)
	{
		end_pos = start_pos + ref_length - 1;
	}
	// insertions and SNPs do not have length >1

	interval_t interval = boost::icl::discrete_interval<pos_t>::right_open(start_pos, end_pos);
	value_t record_set;
	record_set.insert(vcfRecord.line);
	
	// std::cout << "Adding record {" << interval << ", [" << vcfRecord.line << "]}\n";

	_records += std::make_pair(interval, record_set);
}

void 
RecordTracker::
dump(
    std::ostream& os) const
{
	os << "RecordTracker\n";
	for (record_t::const_iterator it = _records.begin(); it != _records.end(); ++it)
	{
		os << it->first << " => {";
		for (value_t::const_iterator vit = it->second.begin(); vit != it->second.end(); ++vit)
		{
			os << *vit << ",";
		}
		os << "}\n";
	}
}

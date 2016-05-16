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
    std::string& record) const
{
	if (_records.empty()) return false;

	const auto posIter(_records.upper_bound(known_pos_range2(beginPos, beginPos)));
	if (posIter == _records.end()) return false;

	if ((*posIter).first.begin_pos() < endPos)
	{
		record = (*posIter).second;
		return true;
	}
	return false;
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

	const known_pos_range2 record_key(start_pos, end_pos);
	
	std::pair<std::map<known_pos_range2, std::string>::iterator, bool> ret;
	ret = _records.insert(std::pair<const known_pos_range2, std::string>(record_key, vcfRecord.line));
	assert(ret.second); // assuming no duplicate records are possible
}

void 
RecordTracker::
dump(
    std::ostream& os) const
{
	os << "RecordTracker\n";
	for (const auto& val : _records)
	{
		os << val.first << " => " << val.second << std::endl;
	}
}

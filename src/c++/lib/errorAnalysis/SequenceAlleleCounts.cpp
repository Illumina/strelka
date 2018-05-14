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

#include "SequenceAlleleCounts.hh"
#include "common/Exceptions.hh"
#include "boost/archive/binary_iarchive.hpp"
#include "boost/archive/binary_oarchive.hpp"

#include <fstream>
#include <iostream>



void
SequenceAlleleCounts::
merge(
    const SequenceAlleleCounts& in)
{
    if (!_sampleName.empty() && !in._sampleName.empty())
    {
        if (_sampleName != in._sampleName)
        {
            using namespace illumina::common;
            std::ostringstream oss;
            oss << "Attempted to merge SequenceAlleleCounts with different sample names: '" << _sampleName << "' and '" << in._sampleName << "'";
            BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
        }
    }
    _sampleName = in._sampleName;
    _bases.merge(in._bases);
    _indels.merge(in._indels);
}


void
SequenceAlleleCounts::
save(
    const char* filename) const
{
    using namespace boost::archive;

    assert(nullptr != filename);
    std::ofstream ofs(filename, std::ios::binary);
    binary_oarchive oa(ofs);

    oa << _sampleName;
    oa << _bases;
    oa << _indels;
}



void
SequenceAlleleCounts::
load(
    const char* filename)
{
    using namespace boost::archive;

    clear();

    assert(nullptr != filename);
    std::ifstream ifs(filename, std::ios::binary);
    binary_iarchive ia(ifs);

    ia >> _sampleName;
    ia >> _bases;
    ia >> _indels;
}



void
SequenceAlleleCounts::
dump(
    std::ostream& os) const
{
    getBaseCounts().dump(os);
    getIndelCounts().dump(os);
}

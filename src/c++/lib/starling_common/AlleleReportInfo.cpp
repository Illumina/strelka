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

#include "starling_common/AlleleReportInfo.hh"

#include <iostream>
#if 0
#include "blt_common/ref_context.hh"
#include "blt_util/seq_util.hh"
#include "starling_common/pos_basecall_buffer.hh"

#include "boost/lexical_cast.hpp"

#include <cassert>
#include <cmath>
#endif



void AlleleSampleReportInfo::dump(std::ostream& os) const
{
    os << "n_confident_ref_reads=" << n_confident_ref_reads
       << ",n_confident_indel_reads=" << n_confident_indel_reads
       << ",n_confident_alt_reads=" << n_confident_alt_reads
       << ",n_other_reads=" << n_other_reads
       << ",indelLocusDepth=" << indelLocusDepth;
}

std::ostream& operator<<(std::ostream& os, const AlleleSampleReportInfo& obj)
{
    obj.dump(os);
    return os;
}

void AlleleReportInfo::dump(std::ostream& os) const
{
    os << "repeatUnit=" << repeatUnit
       << ",ref_repeat_count=" << refRepeatCount
       << ",indel_repeat_count=" << indelRepeatCount
       << ",interruptedHomopolymerLength=" << interruptedHomopolymerLength
       << ",it=" << it;
}

std::ostream& operator<<(std::ostream& os, const AlleleReportInfo& obj)
{
    obj.dump(os);
    return os;
}

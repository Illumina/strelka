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


#include "bam_seq_read_util.hh"



unsigned
getReadAmbiguousEndLength(
    const bam_seq& bseq,
    const bool isFwdStrand)
{
    if (isFwdStrand)
    {
        unsigned read_end(bseq.size());
        while ((read_end>0) && (bseq.get_char(read_end-1)=='N'))
        {
            read_end--;
        }
        return (bseq.size()-read_end);
    }
    else
    {
        const unsigned bsize(bseq.size());
        unsigned read_start(0);
        while ((read_start<bsize) && (bseq.get_char(read_start)=='N'))
        {
            read_start++;
        }
        return read_start;
    }
}

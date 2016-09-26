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

/// \file

/// \author Chris Saunders
///

#include "blt_util/log.hh"
#include "blt_util/seq_util.hh"
#include "starling_common/starling_ref_seq.hh"

extern "C" {
#include "htslib/faidx.h"
}

#include <cassert>
#include <cstdlib>

#include <iostream>
#include <sstream>


static
void
get_samtools_ref_seq(const char* ref_file,
                     const char* chr_name,
                     std::string& ref_seq,
                     const pos_range& pr)
{
    faidx_t* fai(fai_load(ref_file));
    std::ostringstream fa_region_oss;
    fa_region_oss << chr_name;
    if (pr.is_end_pos)
    {
        const pos_t begin(1+(pr.is_begin_pos ? pr.begin_pos : 0));
        fa_region_oss << ':' << begin << '-' << pr.end_pos;
    }
    else if (pr.is_begin_pos)
    {
        fa_region_oss << ':' << pr.begin_pos+1;
    }
    int len; // throwaway...
    char* ref_tmp(fai_fetch(fai,fa_region_oss.str().c_str(), &len));
    if (NULL == ref_tmp)
    {
        log_os << "ERROR: Can't find sequence region '" << fa_region_oss.str() << "' in reference file: '" << ref_file << "'\n";
        exit(EXIT_FAILURE);
    }
    ref_seq.assign(ref_tmp);
    free(ref_tmp);
    fai_destroy(fai);
}



void
get_starling_ref_seq(const starling_base_options& opt,
                     reference_contig_segment& ref)
{
    assert(not opt.referenceFilename.empty());

    // make a temp copy of report_range here to determine how much we
    // pull from ref_seq:
    pos_range ref_range = opt.user_report_range;
    static const pos_t region_read_size_pad(512);
    const pos_t pad_size(opt.max_indel_size+region_read_size_pad);
    if (ref_range.is_begin_pos)
    {
        ref_range.begin_pos -= 1;
        ref_range.begin_pos = std::max(0,ref_range.begin_pos-pad_size);
    }
    else
    {
        ref_range.set_begin_pos(0);
    }

    if (ref_range.is_end_pos)
    {
        ref_range.end_pos += pad_size;
    }

    ref.set_offset(ref_range.begin_pos);

    assert(! opt.bam_seq_name.empty());
    get_samtools_ref_seq(opt.referenceFilename.c_str(), opt.bam_seq_name.c_str(), ref.seq(), ref_range);
    standardize_ref_seq(opt.referenceFilename.c_str(), opt.bam_seq_name.c_str(), ref.seq(), ref_range.begin_pos);
}

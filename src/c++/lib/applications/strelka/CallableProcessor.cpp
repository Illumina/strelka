// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

///
/// \author Chris Saunders
///

#include "strelka/CallableProcessor.hh"

#include <iostream>



void
CallableProcessor::
add(
    const std::string& chrom,
    const pos_t outputPos,
    const somatic_snv_genotype_grid& sgtg)
{
    if ((sgtg.rs.snv_qphred < _minQSS) && (sgtg.rs.nonsomatic_qphred < _minNQSS)) return;

    // determine if we need to flush current range:
    if (_is_range)
    {
        if ((chrom != _chrom) || ((_prange.end_pos+1) != outputPos))
        {
            flush();
        }
    }

    if (_is_range)
    {
        _prange.set_end_pos(outputPos);
    }
    else
    {
        _chrom=chrom;
        _prange.set_begin_pos(outputPos-1);
        _prange.set_end_pos(outputPos);
        _is_range=true;
    }
}



void
CallableProcessor::
flush()
{
    if (NULL == _osptr) return;
    if (! _is_range) return;

    (*_osptr) << _chrom << '\t' << _prange.begin_pos << '\t' << _prange.end_pos << '\n';

    _is_range=false;
}

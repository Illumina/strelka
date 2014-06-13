// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Starka
// Copyright (c) 2009-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

/// \file

/// \author Chris Saunders
///

#pragma once

#include "blt_util/bam_seq.hh"
#include "blt_util/seq_util.hh"

#include "boost/array.hpp"

#include "stdint.h"

#include <iosfwd>
#include <vector>


/// the next two structures are used to generate the GATK haplotype score:
///
struct hap_cand
{

    hap_cand(const bam_seq_base& read_seq,
             const uint8_t* init_qual,
             const int offset);  // the offset into read of the pileup base

    uint16_t
    total_qual() const
    {
        return _total_qual;
    }

    uint8_t
    base_id(const unsigned i) const
    {
        assert(i<HAP_SIZE);
        return (_bq[i] ? (_bq[i] & BASE_MASK) : BASE_ID::ANY );
    }

    uint8_t
    qual(const unsigned i) const
    {
        assert(i<HAP_SIZE);
        return (_bq[i]>>QUAL_SHIFT);
    }

    bool
    operator<(const hap_cand& rhs) const
    {
        // sort the highest scores first:
        return (_total_qual > rhs._total_qual);
    }

    enum { FLANK_SIZE = 10,
           HAP_SIZE = FLANK_SIZE*2+1,
           BASE_MASK = 0x3,
           QUAL_SHIFT = 2
         };

private:
    uint16_t _total_qual;
    boost::array<uint8_t,HAP_SIZE> _bq;
};


std::ostream& operator<<(std::ostream& os, const hap_cand& hc);


typedef std::vector<hap_cand> hap_set_t;


double
get_hapscore(hap_set_t& hap_set);

// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/downloads/sequencing/licenses/>.
//

/// \file

/// \author Chris Saunders
///
#ifndef __VCF_RECORD_HH
#define __VCF_RECORD_HH

#include <iosfwd>
#include <string>
#include <vector>


struct vcf_record {

    vcf_record() : pos(0) { clear(); }

    // set record from record string s, return false on error
    bool set(const char* s);

    void clear() {
        chrom="";
        pos=0;
        ref="";
        alt.clear();
    }

    std::string chrom;
    int pos;
    std::string ref;
    std::vector<std::string> alt;
};


std::ostream& operator<<(std::ostream& os, const vcf_record& vcfr);


#endif

// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2012 Illumina, Inc.
//
// This software is covered by the "Illumina Genome Analyzer Software
// License Agreement" and the "Illumina Source Code License Agreement",
// and certain third party copyright/licenses, and any user of this
// source file is bound by the terms therein (see accompanying files
// Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
// Illumina_Source_Code_License_Agreement.pdf and third party
// copyright/license notices).
//
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
    bool set(const char* s,
             const int len);

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

private:
    enum { MAX_WORD_COUNT = 6 };
    char* _word[MAX_WORD_COUNT];
};


std::ostream& operator<<(std::ostream& os, const vcf_record& vcfr);


#endif

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

#include "blt_util/parse_util.hh"
#include "blt_util/vcf_record.hh"

#include <cassert>
#include <cctype>

#include <algorithm>
#include <iostream>



struct convert {
    void operator()(char& c) { c = toupper((unsigned char)c); }
};


static
void
stoupper(std::string& s) {
    std::for_each(s.begin(), s.end(), convert());
}



bool
vcf_record::
set(const char* s,
    const int len) {

    static const char sep('\t');
    static const unsigned maxword(5);

    // simple tab parse:
    const char* start(s);
    const char* p(start);

    unsigned wordindex(0);
    while(wordindex<maxword) {
        if ((*p == sep) || (*p == '\n') || (*p == '\0')) {
            switch(wordindex) {
            case 0:
                chrom=std::string(start,p-start);
                break;
            case 1:
                pos=casava::blt_util::parse_int(start);
                assert(start==p);
                break;
            case 2:
                // skip this field...
                break;
            case 3:
                ref=std::string(start,p-start);
                stoupper(ref);
                break;
            case 4:
                // addition parse loop for ',' character:
                {
                    const char* p2(start);
                    while(p2<=p){
                        if((*p2==',') || (p2==p)) {
                            alt.push_back(std::string(start,p2-start));
                            stoupper(alt.back());
                            start=p2+1;
                        }
                        p2++;
                    }
                }
                break;
            default:
                assert(0);
            }
            start=p+1;
            wordindex++;
        }
        if((*p == '\n') || (*p == '\0')) break;
        ++p;
    }

    return (wordindex >= maxword);
}



std::ostream& operator<<(std::ostream& os, const vcf_record& vcfr) {

    os << vcfr.chrom << '\t'
       << vcfr.pos << '\t'
       << '.' << '\t'
       << vcfr.ref << '\t';

    const unsigned nalt(vcfr.alt.size());
    for(unsigned a(0);a<nalt;++a) {
        if(a) os << ',';
        os << vcfr.alt[a];
    }
    os << '\t'
       << '.' << '\t'
       << '.' << '\t'
       << '.' << '\t'
       << '.' << '\n';

    return os;
}


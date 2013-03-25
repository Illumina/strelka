// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2012 Illumina, Inc.
//
// Permission is hereby granted, free of charge, to any person
// obtaining a copy of this software and associated documentation
// files (the "Software"), to deal in the Software without
// restriction, including without limitation the rights to use, copy,
// modify, merge, publish, distribute, sublicense, and/or sell copies
// of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
// BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
// ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
// CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
//

/// \file
/// random vcf utilities

/// \author Chris Saunders
///
#ifndef __VCF_UTIL_HH
#define __VCF_UTIL_HH


#include <cstring>
#include <iosfwd>
#include <vector>


namespace VCFID {
enum index_t
{
    CHROM,
    POS,
    ID,
    REF,
    ALT,
    QUAL,
    FILT,
    INFO,
    FORMAT,
    SAMPLE,
    SIZE
};
}



inline
const char*
vcf_col_label() {
    static const char h[] = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
    return h;
}


std::ostream&
vcf_fileDate(std::ostream& os);



// look for 'key' in vcf FORMAT field, provide index of key or return
// false
//
inline
bool
get_format_key_index(const char* format,
                     const char* key,
                     unsigned& index) {
    index=0;
    do {
        if(index) format++;
        if(0==strncmp(format,key,strlen(key))) return true;
        index++;
    } while(NULL != (format=strchr(format,':')));
    return false;
}



// return pointer to
//
inline
const char*
get_format_string_nocopy(const char* const* word,
                         const char* key) {

    unsigned keynum(0);
    if(! get_format_key_index(word[VCFID::FORMAT],key,keynum)) return NULL;

    const char* sample(word[VCFID::SAMPLE]);
    for(; keynum; sample++) {
        if(! *sample) return NULL;
        if((*sample)==':') keynum--;
    }
    return sample;
}



// returns -1 for '.' alleles
void
parse_gt(const char* gt,
         std::vector<int>& gti,
         const bool is_allow_bad_end_char=false);


#endif

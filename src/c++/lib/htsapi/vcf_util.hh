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

/// random vcf utilities
///
/// \author Chris Saunders
///

#pragma once

#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstring>

#include <iosfwd>
#include <vector>


namespace VCFID
{
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
vcf_col_label()
{
    static const char h[] = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
    return h;
}


std::ostream&
vcf_fileDate(std::ostream& os);


void
write_vcf_filter(
    std::ostream& os,
    const char* id,
    const char* desc);


// look for 'key' in vcf FORMAT field, provide index of key or return
// false
//
inline
bool
get_format_key_index(const char* format,
                     const char* key,
                     unsigned& index)
{
    index=0;
    do
    {
        if (index) format++;
        if (0==strncmp(format,key,strlen(key))) return true;
        index++;
    }
    while (NULL != (format=strchr(format,':')));
    return false;
}



// return pointer to
//
inline
const char*
get_format_string_nocopy(const char* const* word,
                         const char* key)
{
    unsigned keynum(0);
    if (! get_format_key_index(word[VCFID::FORMAT],key,keynum)) return NULL;

    const char* sample(word[VCFID::SAMPLE]);
    for (; keynum; sample++)
    {
        if (! *sample) return NULL;
        if ((*sample)==':') keynum--;
    }
    return sample;
}



// returns -1 for '.' alleles
void
parse_gt(
    const char* gt,
     std::vector<int>& gti,
     const bool is_allow_bad_end_char=false);


struct VcfGenotypeUtil
{
    static
    unsigned
    getGenotypeCount(
        const unsigned ploidy,
        const uint8_t alleleCount)
    {
        if (ploidy == 1)
        {
            return alleleCount;
        }
        else if (ploidy == 2)
        {
            return getGenotypeIndex(0, alleleCount);
        }
        else
        {
            assert(false and "unexpected ploidy");
        }
    }

    /// allele to genotype functions
    ///
    /// ploidy is implied by the number of arguments
    static
    unsigned
    getGenotypeIndex(
        const uint8_t allele0Index)
    {
        return allele0Index;
    }


    static
    unsigned
    getGenotypeIndex(
        const uint8_t allele0Index,
        const uint8_t allele1Index)
    {
        assert(allele0Index<=allele1Index);
        return allele0Index+(allele1Index*(allele1Index+1)/2);
    }

    /// genotype to allele functions
    ///
    /// ploidy is implied by the number of arguments
    static
    void
    getAlleleIndices(
        const unsigned genotypeIndex,
        uint8_t& allele0Index)
    {
        allele0Index = genotypeIndex;
    }

    static
    void
    getAlleleIndices(
        const unsigned genotypeIndex,
        uint8_t& allele0Index,
        uint8_t& allele1Index)
    {
        // from quadratic inversion of genotype index function
        // 1. Function to invert:
        //     g = (a1^2)/2 + a1/2 + a0
        // 2. Solve for a1 >= 0 when a0 == 0
        //     g = (a1^2)/2 + a1/2 + 0
        //     0 = a1^2 + a1 - 2g
        //     a1 = (-1+sqrt(1+4*2G) )/ 2
        //
        allele1Index = std::floor((std::sqrt(1.0+8.0*genotypeIndex) - 1.0) / 2.0);
        allele0Index = genotypeIndex - allele1Index*(allele1Index+1)/2;
        assert(allele0Index<=allele1Index);
    }

    static
    void
    writeGenotype(
        const unsigned ploidy,
        const unsigned genotypeIndex,
        std::ostream& os)
    {
        if (ploidy == 1)
        {
            uint8_t allele0Index;
            getAlleleIndices(genotypeIndex, allele0Index);
            writeGenotype(allele0Index, os);
        }
        else if (ploidy == 2)
        {
            uint8_t allele0Index;
            uint8_t allele1Index;
            getAlleleIndices(genotypeIndex, allele0Index, allele1Index);
            writeGenotype(allele0Index, allele1Index, os);
        }
        else
        {
            assert(false and "Unexpected ploidy value");
        }
    }

    static
    void
    writeGenotype(
        const uint8_t allele0Index,
        std::ostream& os);

    static
    void
    writeGenotype(
        const uint8_t allele0Index,
        const uint8_t allele1Index,
        std::ostream& os);
};

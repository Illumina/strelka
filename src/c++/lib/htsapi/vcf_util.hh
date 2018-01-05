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
///
/// \brief VCF utilities
///
/// \author Chris Saunders


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


/// look for 'key' in vcf FORMAT field, provide index of key or return
/// false
///
inline
bool
get_format_key_index(const char* format,
                     const char* key,
                     unsigned& index)
{
    index=0;
    do
    {
        if (index>0) format++;
        if (0==strncmp(format,key,strlen(key))) return true;
        index++;
    }
    while (nullptr != (format=strchr(format,':')));
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
    if (! get_format_key_index(word[VCFID::FORMAT],key,keynum)) return nullptr;

    const char* sample(word[VCFID::SAMPLE]);
    for (; keynum>0; sample++)
    {
        if (! *sample) return nullptr;
        if ((*sample)==':') keynum--;
    }
    return sample;
}



/// returns -1 for '.' alleles
void
parse_gt(
    const char* gt,
    std::vector<int>& gti,
    const bool is_allow_bad_end_char=false);


/// a general-purpose genotype object, capable of expressing any valid GT in the VCF standard (of ploidy 2 or less)
struct VcfGenotype
{
    VcfGenotype()
    {
        setGenotypeFromAlleleIndices();
    }

    explicit
    VcfGenotype(const unsigned allele0Index)
    {
        setGenotypeFromAlleleIndices(allele0Index);
    }

    VcfGenotype(
        const unsigned allele0Index,
        const unsigned allele1Index,
        const bool isPhased = false)
    {
        setGenotypeFromAlleleIndices(allele0Index, allele1Index, isPhased);
    }

    void
    setGenotypeFromAlleleIndices()
    {
        _ploidy=0;
    }

    void
    setGenotypeFromAlleleIndices(
        const unsigned allele0Index)
    {
        _ploidy=1;
        _allele0Index = allele0Index;
    }

    void
    setGenotypeFromAlleleIndices(
        const unsigned allele0Index,
        const unsigned allele1Index,
        const bool isPhased = false)
    {
        _ploidy=2;
        _allele0Index = allele0Index;
        _allele1Index = allele1Index;
        _isPhased = isPhased;

        if (not _isPhased)
        {
            if (_allele0Index > _allele1Index)
            {
                std::swap(_allele0Index,_allele1Index);
            }
        }
    }

    void
    setAllele0HaplotypeId(uint8_t haplotypeId)
    {
        _allele0HaplotypeId = haplotypeId;
    }

    void
    setAllele1HaplotypeId(uint8_t complexAlleleId)
    {
        _allele1HaplotypeId = complexAlleleId;
    }

    void
    addAltAlleleHaplotypeCountRatio(const float altAlleleHaplotypeCountRatio)
    {
        _altAlleleHaplotypeCountRatio += altAlleleHaplotypeCountRatio;
    }

    void
    setPhased(
        const bool isFlip = false)
    {
        _isPhased = true;
        if (not isFlip) return;
        std::swap(_allele0Index,_allele1Index);
    }

    bool
    isUnknown() const
    {
        return ((_ploidy < 1) || (_ploidy > 2));
    }

    bool
    getIsPhased() const
    {
        return _isPhased;
    }

    bool
    isHet() const
    {
        if (getPloidy() != 2)
        {
            return false;
        }
        return (getAllele0Index() != getAllele1Index());
    }

    bool
    containsReference() const
    {
        bool isReferenceContained = (getAllele0Index() == 0);
        if (!isReferenceContained and (getPloidy() == 2))
        {
            isReferenceContained = (getAllele1Index() == 0);
        }
        return isReferenceContained;
    }

    uint8_t
    getAllele0Index() const
    {
        assert(_ploidy >= 1);
        return _allele0Index;
    }

    uint8_t
    getAllele1Index() const
    {
        assert(_ploidy >= 2);
        return _allele1Index;
    }

    uint8_t
    getAllele0HaplotypeId() const
    {
        assert(_ploidy >= 1);
        return _allele0HaplotypeId;
    }

    uint8_t
    getAllele1HaplotypeId() const
    {
        assert(_ploidy >= 2);
        return _allele1HaplotypeId;
    }

    float
    getAltHaplotypeCountRatio() const
    {
        // _altAlleleHaplotypeCountRatio may exceed 1.0
        // because some reads can be double counted if haplotypes are reconstructed by assembly
        return _altAlleleHaplotypeCountRatio > 1.0f ? 1.0f : _altAlleleHaplotypeCountRatio;
    }

    bool isConflict() const
    {
        // haplotypeId==3 should be hom
        if (isHet())
            return (_allele0HaplotypeId == _allele1HaplotypeId) or (_allele0HaplotypeId == 3) or (_allele1HaplotypeId == 3);

        return _allele0HaplotypeId != _allele1HaplotypeId;
    }

    int
    getPloidy() const
    {
        return _ploidy;
    }

    bool
    isVariant() const
    {
        return ((_allele0Index != 0) || (_allele1Index != 0));
    }

    void
    clear()
    {
        _ploidy = 2;
        _allele0Index = 0;
        _allele1Index = 0;

        _allele0HaplotypeId = 0;
        _allele1HaplotypeId = 0;

        _isPhased = false;
    }

    bool
    operator==(const VcfGenotype& rhs) const
    {
        return
            ((_ploidy == rhs._ploidy) and
             (_allele0Index == rhs._allele0Index) and
             (_allele1Index == rhs._allele1Index) and
             (_isPhased == rhs._isPhased));
    }

private:
    /// these first four values are private because they must change in sync:
    int _ploidy = 2;
    uint8_t _allele0Index = 0;
    uint8_t _allele1Index = 0;

    uint8_t _allele0HaplotypeId = 0;
    uint8_t _allele1HaplotypeId = 0;

    float _altAlleleHaplotypeCountRatio = 0.0f;

    bool _isPhased = false;
};



std::ostream&
operator<<(std::ostream& os, const VcfGenotype& vcfGt);


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
            return 0;
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
        const int ploidy,
        const unsigned genotypeIndex,
        std::ostream& os)
    {
        if (ploidy == 1)
        {
            uint8_t allele0Index;
            getAlleleIndices(genotypeIndex, allele0Index);
            writeGenotype(allele0Index, os);
        }
        else if ((ploidy == 2) or (ploidy == -1))
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

    static
    void
    writeGenotype(
        const VcfGenotype& vcfGt,
        std::ostream& os);
};


inline
void
setVcfGenotypeFromGenotypeIndex(
    const int genotypePloidy,
    const unsigned genotypeIndex,
    VcfGenotype& vcfGt)
{
    if (genotypePloidy == 1)
    {
        uint8_t allele0Index;
        VcfGenotypeUtil::getAlleleIndices(genotypeIndex, allele0Index);
        vcfGt.setGenotypeFromAlleleIndices(allele0Index);
    }
    else if (genotypePloidy == 2)
    {
        uint8_t allele0Index, allele1Index;
        VcfGenotypeUtil::getAlleleIndices(genotypeIndex, allele0Index, allele1Index);
        vcfGt.setGenotypeFromAlleleIndices(allele0Index, allele1Index);
    }
    else
    {
        assert(false and "Unexpected ploidy");
    }
}

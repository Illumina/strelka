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


#include "blt_common/ref_context.hh"
#include "blt_util/seq_util.hh"
#include <cassert>



unsigned
get_left_shifted_hpol_size(
    const pos_t pos,
    const reference_contig_segment& ref)
{
    unsigned size(1);

    const char rpos(ref.get_base(pos));
    if (pos>0)
    {
        if (rpos == ref.get_base(pos-1)) return size;
    }

    const pos_t rs(ref.end());
    for (pos_t i(pos+1); i<rs; i++)
    {
        if (rpos != ref.get_base(i)) break;
        size++;
    }
    return size;
}

bool
isLeftEndOfSTR(
    const unsigned patternSize,
    const pos_t pos,
    const reference_contig_segment& ref)
{
    // the current position cannot be the start of an STR track
    // if the previous base is the same as the end of the potential pattern
    return ref.get_base(pos - 1) != ref.get_base(pos + patternSize - 1);
}

unsigned
getLeftShiftedSTRRepeatCount(
    const unsigned patternSize,
    const pos_t startPosition,
    const reference_contig_segment& ref)
{
    unsigned size(1);

    if (startPosition > (int) patternSize - 1)
    {
        if (!isLeftEndOfSTR(patternSize, startPosition, ref))
        {
            return size;
        }
    }

    const pos_t endOfRef(ref.end());
    for (pos_t posInRef(startPosition + patternSize); posInRef < endOfRef; posInRef += patternSize)
    {
        if (!compareRepeatPattern(patternSize, startPosition, posInRef, ref))
        {
            break;
        }
        size++;
    }
    return size;
}

bool
compareRepeatPattern(
    const unsigned repeatPatternSize,
    const unsigned pos1,
    const unsigned pos2,
    const reference_contig_segment& ref)
{
    for (int posInPattern(repeatPatternSize -1 ); posInPattern >= 0; posInPattern--)
    {
        if (ref.get_base(pos1+posInPattern) != ref.get_base(pos2 + posInPattern))
        {
            return false;
        }
    }
    return true;
}

void searchForSTR(
    const unsigned patternSize,
    const unsigned pos,
    bool& isBaseInStr,
    bool& isBaseLeftEndOfStr,
    const reference_contig_segment& ref)
{
    isBaseInStr = false;
    isBaseLeftEndOfStr = false;

    // start from the k-mer at pos vs. k-mer at pos+k
    // end at k-mer at pos-k-patternSize+1 vs. k-mer at pos-patternSize+1
    for (unsigned p = 0; p < patternSize; p++)
    {
        if (compareRepeatPattern(patternSize, pos - p, pos + patternSize - p, ref) ||
            compareRepeatPattern(patternSize, pos - p, pos - patternSize - p, ref))
        {
            if (0 == p && isLeftEndOfSTR(patternSize, pos, ref))
            {
                isBaseLeftEndOfStr = true;
            }
            isBaseInStr = true;
            return;
        }
    }

}

unsigned
get_snp_hpol_size(const pos_t pos,
                  const reference_contig_segment& ref)
{
    // count upstream repeats:
    bool is_up_repeat(false);
    char up_repeat('N');
    unsigned up_size(0);
    for (pos_t i(pos-1); i>=0; i--)
    {
        if (is_up_repeat)
        {
            if (up_repeat != ref.get_base(i)) break;
        }
        else
        {
            up_repeat=ref.get_base(i);
            is_up_repeat=true;
            if (up_repeat=='N') break;
        }
        up_size++;
    }

    // count downstream repeats:
    bool is_dn_repeat(false);
    char dn_repeat('N');
    unsigned dn_size(0);
    const pos_t rs(ref.end());
    for (pos_t i(pos+1); i<rs; i++)
    {
        if (is_dn_repeat)
        {
            if (dn_repeat != ref.get_base(i)) break;
        }
        else
        {
            dn_repeat=ref.get_base(i);
            is_dn_repeat=true;
            if (dn_repeat=='N') break;
        }
        dn_size++;
    }

    return 1+((dn_repeat==up_repeat) ? up_size+dn_size : std::max(up_size,dn_size) );
}



// helper class for finding interrupted homopolymers:
//
struct ihpol_data
{
    ihpol_data()
        : r1('N'), r2('N')
        , nr1(0), nr2(0)
        , size(0)
    {}

    // returns false when limit is reached:
    bool
    add_base(const char b)
    {
        if       (nr1==0)
        {
            r1=b;
            nr1=1;
        }
        else if (r1==b)
        {
            if ((nr2>1) || (r1=='N')) return false;
            nr1++;
        }
        else if (nr2==0)
        {
            r2=b;
            nr2=1;
        }
        else if (r2==b)
        {
            if ((nr1>1) || (r2=='N')) return false;
            nr2++;
        }
        else
        {
            return false;
        }
        size++;
        return true;
    }

    unsigned
    max_allele_size() const
    {
        return std::max(nr1,nr2);
    }

    char r1,r2;
    unsigned nr1,nr2;
    unsigned size;
};



unsigned
getInterruptedHomopolymerLength(
    const pos_t pos,
    const reference_contig_segment& ref)
{
    // count current base + upstream repeats:
    ihpol_data up_ihd;
    for (pos_t i(pos); i>=0; i--)
    {
        const char base(ref.get_base(i));
        if (! up_ihd.add_base(base)) break;
    }
    // then see how far we can extend downstream:
    const pos_t rs(ref.end());
    for (pos_t i(pos+1); i<rs; i++)
    {
        const char base(ref.get_base(i));
        if (! up_ihd.add_base(base)) break;
    }

    // count current base + downstream repeats:
    ihpol_data dn_ihd;
    for (pos_t i(pos); i<rs; i++)
    {
        const char base(ref.get_base(i));
        if (! dn_ihd.add_base(base)) break;
    }
    // then see how far we can extend upstream:
    for (pos_t i(pos-1); i>=0; i--)
    {
        const char base(ref.get_base(i));
        if (! dn_ihd.add_base(base)) break;
    }

    // return the longer of the two:
    return std::max(up_ihd.max_allele_size(),dn_ihd.max_allele_size());
}

unsigned
shortestUnencountered(
    const reference_contig_segment& ref,
    const pos_t pos,
    const unsigned numEncoded,
    const bool left)
{
    std::string encoded;
    std::string newPrefix;
    // get substring encoded thus far:
    if (left)
    {
        ref.get_substring(pos-numEncoded+1,numEncoded,encoded);
    }
    else
    {
        ref.get_substring(pos,numEncoded,encoded);
    }
    unsigned len(0);
    do
    {
        ++len;
        if (left)
        {
            ref.get_substring(pos-numEncoded-len+1,len,newPrefix);
        }
        else
        {
            ref.get_substring(pos+numEncoded,len,newPrefix);
        }
    }
    // keep looking while string of size len has already been encountered:
    while (encoded.find(newPrefix) != std::string::npos);
    return len;
}

unsigned
computeContextCompressability(
    const reference_contig_segment& ref,
    const pos_t leftPos,
    const pos_t rightPos,
    const unsigned numKeys)
{
    assert(numKeys>0);
    assert(leftPos <= rightPos);

    // complexity to the left:
    unsigned numEncodedLeft(1);
    unsigned numEncodedRight(1);
    for (unsigned key(0); key<numKeys-1; ++key)
    {
        // complexity to the left:
        numEncodedLeft += shortestUnencountered(ref,leftPos-1,numEncodedLeft,true);
        // complexity to the right:
        numEncodedRight += shortestUnencountered(ref,rightPos,numEncodedRight,false);
    }

    return std::max(numEncodedLeft, numEncodedRight);
}

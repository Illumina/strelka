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

///
/// \author Chris Saunders
///

#include "vcf_record.hh"
#include "blt_util/parse_util.hh"

#include <cassert>
#include <cctype>

#include <algorithm>
#include <iostream>



struct convert
{
    void operator()(char& c) const
    {
        c = toupper((unsigned char)c);
    }
};


static
void
stoupper(std::string& s)
{
    std::for_each(s.begin(), s.end(), convert());
}



bool
vcf_record::
set(const char* s)
{
    static const char sep('\t');
    static const unsigned maxword(5);

    clear();

    line = s;

    // simple tab parse:
    const char* start(s);
    const char* p(start);

    unsigned wordindex(0);
    while (wordindex<maxword)
    {
        if ((*p == sep) || (*p == '\n') || (*p == '\0'))
        {
            switch (wordindex)
            {
            case 0:
                chrom=std::string(start,p-start);
                break;
            case 1:
                pos=illumina::blt_util::parse_int(start);
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
                // additional parse loop for ',' character:
            {
                const char* p2(start);
                while (p2<=p)
                {
                    if ((*p2==',') || (p2==p))
                    {
                        alt.emplace_back(start,p2-start);
                        stoupper(alt.back());
                        start=p2+1;
                    }
                    p2++;
                }
            }
            break;
            default:
                assert(0);
                break;
            }
            start=p+1;
            wordindex++;
        }
        if ((*p == '\n') || (*p == '\0')) break;
        ++p;
    }

    return (wordindex >= maxword);
}

bool
vcf_record::
is_normalized() const
{
    // normalized indels are left-aligned, reference-padded, and parsimonious
    // normalized SNVs are a single differing base
    // normalized MNVs (and complex alleles) have differing bases at the beginning
    // and end of the ref and alt alleles
    // see http://genome.sph.umich.edu/wiki/Variant_Normalization
    unsigned ref_length = ref.size();
    assert (ref_length != 0);

    for (const auto& alt_allele : alt)
    {
        // all normalized variants must differ at the last ref and alt base
        // this checks for right-padding, i.e. parsimony
        if((*alt_allele.rbegin()) == (*ref.rbegin()))
        {
            return false;
        }

        unsigned alt_length = alt_allele.size();
        assert (alt_length != 0);
        if(alt_length != ref_length)
        {
            // this checks that indels are reference-padded
            if ( (*alt_allele.begin()) == (*ref.begin()))
            {
                // this checks that they're left-shifted
                for (unsigned i = ref_length - 1, j = alt_length - 1; ; ++i, ++j)
                {
                    if (ref[i] != alt_allele[j])
                    {
                        break;
                    }
                    else if (i == 0 || j == 0)
                    {
                        return false;
                    }
                }
            }
            // if the first and last bases of alleles with two differing lengths
            // do not match, the record represents a complex allele, and fulfills
            // normalization requirements
        }
        else
        {
            // we already checked that the last base of an SNV/MNV
            // differed between alt and ref at the beginning of the
            // loop
            if (ref_length > 1)
            {
                if ((*alt_allele.begin()) == (*ref.begin()))
                {
                    return false;
                }
            }
        }
    }
    return true;
}

std::ostream& operator<<(std::ostream& os, const vcf_record& vcfr)
{
    os << vcfr.chrom << '\t'
       << vcfr.pos << '\t'
       << '.' << '\t'
       << vcfr.ref << '\t';

    const unsigned nalt(vcfr.alt.size());
    for (unsigned a(0); a<nalt; ++a)
    {
        if (a) os << ',';
        os << vcfr.alt[a];
    }
    os << '\t'
       << '.' << '\t'
       << '.' << '\t'
       << '.' << '\t'
       << '.' << '\n';

    return os;
}


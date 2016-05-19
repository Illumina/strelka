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

#pragma once

#include "blt_util/qscore.hh"

#include "boost/utility.hpp"

#include <iosfwd>
#include <vector>


namespace STAR_DIINDEL
{
enum index_t
{
    NOINDEL,
    HOM,
    HET,
    SIZE
};

inline
const char*
label(const unsigned idx)
{
    switch (idx)
    {
    case NOINDEL:
        return "ref";
    case HOM:
        return "hom";
    case HET:
        return "het";
    default:
        return "xxx";
    }
}

inline
const char*
get_gt_label(const unsigned idx)
{
    switch (idx)
    {
    case NOINDEL:
        return "0/0";
    case HOM:
        return "1/1";
    case HET:
        return "0/1";
    default:
        assert(false && "Unknown Indel GT");
        return nullptr;
    }
}

inline
uint8_t
get_allele(
    const unsigned idx,
    const unsigned chromidx)
{
    assert(idx<SIZE);
    assert(chromidx<2);
    static const uint8_t res[SIZE][2] = {{0,0},{1,1},{0,1}};
    return res[idx][chromidx];
}


#if 0
// states are the number of copies of I,R,NR
// I=called indel allele
// R=reference (no indels which interfere with I)
// NR=first indel allele interfering with I
enum index_t
{
    I2,
    I1R1,
    I1NR1,
    R2,
    R1NR1,
    NR2,
    SIZE
};

inline
const char*
label(const unsigned idx)
{
    switch (idx)
    {
    case I2:
        return "hom";
    case I1R1:
        return "het";
    case I1NR1:
        return "nonref_het";
    case R2:
        return "ref";
    case R1NR1:
    case NR2:
        return "other";
    default:
        return "xxx";
    }
}
#endif
}


// separate out core with automatic copy ctor:
//
struct GermlineDiploidIndelSimpleGenotypeInfoCore
{
    GermlineDiploidIndelSimpleGenotypeInfoCore()
        : is_indel(false), ploidy(2), max_gt(0), max_gt_poly(0), phredLoghood(STAR_DIINDEL::SIZE,0)
    {
        static const int qp(error_prob_to_qphred((1.-init_p())));
        indel_qphred=qp;
        max_gt_qphred=qp;
        max_gt_poly_qphred=qp;
    }

    bool
    is_haploid() const
    {
        return (1 == ploidy);
    }

    bool
    is_noploid() const
    {
        return (0 == ploidy);
    }

    // debug
    void dump(std::ostream& os) const;

    /// only applies to PLs so far:
    static const int maxQ;

protected:

    static
    double
    init_p()
    {
        static const double p(1./static_cast<double>(STAR_DIINDEL::SIZE));
        return p;
    }

public:

    bool is_indel;

    // hack haploid/'noploid' model into diploid data structure:
    // only {2,1,0} are actually used at present
    int ploidy;

    unsigned max_gt;
    int indel_qphred;
    int max_gt_qphred;

    unsigned max_gt_poly;
    int max_gt_poly_qphred;

    bool is_forced_output;
    bool is_zero_coverage;

    std::vector<unsigned> phredLoghood;
};



struct starling_diploid_indel : public GermlineDiploidIndelSimpleGenotypeInfoCore, private boost::noncopyable
{
    starling_diploid_indel()
        : GermlineDiploidIndelSimpleGenotypeInfoCore()
    {
        static const double p(init_p());
        for (unsigned i(0); i<STAR_DIINDEL::SIZE; ++i) pprob[i] = p;
    }

    void dump(std::ostream& os) const;

    // members:
    double pprob[STAR_DIINDEL::SIZE];
};


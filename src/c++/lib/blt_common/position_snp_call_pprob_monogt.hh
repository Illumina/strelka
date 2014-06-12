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
#ifndef __POSITION_SNP_CALL_PPROB_MONOGT_HH
#define __POSITION_SNP_CALL_PPROB_MONOGT_HH

#include "blt_common/snp_pos_info.hh"

#include "blt_util/monogt.hh"

#include <boost/utility.hpp>

#include <iosfwd>



struct monoploid_genotype : private boost::noncopyable {

    monoploid_genotype() : is_snp(false), ref_gt(0), max_gt(0), max2_gt(0) {
        static const double p(1./static_cast<double>(MONOGT::SIZE));
        for (unsigned i(0); i<MONOGT::SIZE; ++i) pprob[i] = p;
    }

    bool is_snp;
    unsigned ref_gt;
    unsigned max_gt;
    unsigned max2_gt;
    double pprob[MONOGT::SIZE];
};

std::ostream& operator<<(std::ostream& os,monoploid_genotype& dp);



/// \brief call a snp @ pos by calculating the posterior probability
/// of all possible genotypes for a monoploid individual.
///
void
position_snp_call_pprob_monogt(const double theta,
                               const snp_pos_info& pi,
                               monoploid_genotype& dgt);

#endif

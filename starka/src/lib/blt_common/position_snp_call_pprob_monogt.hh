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
#ifndef __POSITION_SNP_CALL_PPROB_MONOGT_HH
#define __POSITION_SNP_CALL_PPROB_MONOGT_HH

#include "blt_common/snp_pos_info.hh"

#include "blt_util/monogt.hh"

#include <boost/utility.hpp>

#include <iosfwd>



struct monoploid_genotype : private boost::noncopyable {

    monoploid_genotype() : is_snp(false), ref_gt(0), max_gt(0), max2_gt(0) {
        static const double p(1./static_cast<double>(MONOGT::SIZE));
        for(unsigned i(0);i<MONOGT::SIZE;++i) pprob[i] = p;
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

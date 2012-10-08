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
#ifndef __POSITION_SNP_CALL_PPROB_NPLOID_HH
#define __POSITION_SNP_CALL_PPROB_NPLOID_HH

#include "blt_common/snp_pos_info.hh"

#include "blt_util/nploid_genotype_util.hh"

#include <iosfwd>


struct nploid_genotype {

    nploid_genotype(const nploid_info& ni)
        : is_snp(false), ref_gt(0), max_gt(0), max2_gt(0),
          pprob(ni.gtype_size()), ploidy(ni.ploidy()) {
        const unsigned n_gt(ni.gtype_size());
        static const double p(1./static_cast<double>(n_gt));
        for(unsigned i(0);i<n_gt;++i) pprob[i] = p;
    }

    bool is_snp;
    unsigned ref_gt;
    unsigned max_gt;
    unsigned max2_gt;
    std::vector<double> pprob;
    unsigned ploidy;
};


void
nploid_write(const nploid_info& ninfo,
             const nploid_genotype& ng,
             std::ostream& os);


/// \brief call a snp @ pos by calculating the posterior probability
/// of all possible genotypes for an n-ploid individual sample, ploidy
/// is specified by nploid_info.
///
void
position_snp_call_pprob_nploid(const double snp_prob,
                               const snp_pos_info& pi,
                               const nploid_info& ninfo,
                               nploid_genotype& ngt);

#endif

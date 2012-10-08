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

/// variation on the original strowman snv caller -- implements a
/// compile-time specified grid in allele frequency space and requires
/// similar frequency as definition of non-somatic.
///

/// \author Chris Saunders
///
#ifndef __POSITION_SOMATIC_SNV_GRID_SHARED_HH
#define __POSITION_SOMATIC_SNV_GRID_SHARED_HH


struct somatic_snv_genotype_grid {

    somatic_snv_genotype_grid()
        : is_snv(false) {}

    typedef bool tier_t;

    struct result_set {
        unsigned ntype;
        unsigned max_gt;
        int snv_qphred;
        int snv_from_ntype_qphred;
    };

    bool is_snv;
    tier_t snv_tier;
    tier_t snv_from_ntype_tier;
    unsigned ref_gt;
    result_set rs;
};

#endif

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
///
/// \author Chris Saunders
///

#include "gvcf_aggregator.hh"

#include <iostream>



void
gvcf_aggregator::
add_site(const pos_t output_pos,
         const char ref,
         const unsigned n_used_calls,
         const unsigned n_unused_calls,
         const snp_pos_info& good_pi,
         const diploid_genotype& dgt,
         const bool is_nf_snp,
         const double sb,
         const unsigned hpol) {

}



void
gvcf_aggregator::
add_indel(const pos_t output_pos,
          const starling_diploid_indel& dindel,
          const starling_indel_report_info& iri,
          const starling_indel_sample_report_info& isri) {

}

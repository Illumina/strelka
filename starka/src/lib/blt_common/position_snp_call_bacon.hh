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
#ifndef __POSITION_SNP_CALL_BACON_HH
#define __POSITION_SNP_CALL_BACON_HH

#include "blt_common/snp_pos_info.hh"


/*
Original bacon doc (from perl module):

=head1 NAME

Illumina::Sequencing::BACON4.pm - BAyesian CONsensus caller.

=head1 SYNOPSIS
BAyesian CONsensus caller.

For each base, work out the probability of seeing the observed
coverage pattern assuming the base is NOT observed
Hopefully one will have staggeringly low probability
This is the most probable allele and it gets a log score -
-10 x log10(prob B not present)
Normalised by smallest value found - convert to log odds ratio
If allele with next lowest probability falls below a threshold
then a heterozygote is called

use Illumina::Sequencing::BACON4.pm qw();

=head1 DESCRIPTION

Exports:
    alleleCall

=head1 AUTHORSHIP

Copyright (c) 2006, 2007 Solexa; 2008 Illumina
This source file is covered by the "Illumina Public Source License"
agreement and bound by the terms therein.

Created by A. J. Cox

*/


struct bacon_scores {

    bacon_scores()
        : is_valid(false),
          max_base_id(0), max2_base_id(0),
          max_score(0.), max2_score(0.) {}

    bool is_valid;
    unsigned max_base_id;
    unsigned max2_base_id;
    double max_score;
    double max2_score;
};


/// \brief call a snp @ pos using bacon method
///
void
position_snp_call_bacon(const snp_pos_info& pi,
                        bacon_scores& bas);

#endif

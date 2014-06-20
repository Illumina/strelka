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

///
/// \author Chris Saunders
///

#include "position_somatic_snv.hh"
#include "position_somatic_snv_grid.hh"
#include "position_somatic_snv_strand_grid.hh"
#include "somatic_indel.hh"
#include "somatic_indel_grid.hh"
#include "strelka_shared.hh"



strelka_deriv_options::
strelka_deriv_options(const strelka_options& opt,
                      const reference_contig_segment& ref)
    : base_t(opt,ref)
    , _sscaller_strand_grid(new somatic_snv_caller_strand_grid(opt,pdcaller()))
    , _sicaller_grid(new somatic_indel_caller_grid(opt,incaller()))
{}

/// dtor required to be in the cpp so that unique ptr can access complete data type
strelka_deriv_options::
~strelka_deriv_options() {}

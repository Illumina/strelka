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
/// note coding convention for all ranges '_pos fields' is:
/// XXX_begin_pos is zero-indexed position at the begining of the range
/// XXX_end_pos is zero-index position 1 step after the end of the range
///

#pragma once

#include "strelka_shared.hh"

#include "blt_common/snp_pos_info.hh"

#include <iosfwd>
#include <memory>


struct strelka_pile_caller
{

    strelka_pile_caller(
        strelka_options& opt,
        std::ostream& os);

    void
    call(
        const bool is_somatic_gvcf,
        const unsigned pos,
        snp_pos_info& norm_pi,
        snp_pos_info& tumor_pi);

private:
    strelka_options& _opt;
    std::auto_ptr<strelka_deriv_options> _dopt_ptr;
    std::ostream& _os;
};


void
strelka_pile_test_run(strelka_options& opt);

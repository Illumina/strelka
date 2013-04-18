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
/// note coding convention for all ranges '_pos fields' is:
/// XXX_begin_pos is zero-indexed position at the begining of the range
/// XXX_end_pos is zero-index position 1 step after the end of the range
///


#ifndef STRELKA_PILE_TEST_RUN_
#define STRELKA_PILE_TEST_RUN_

#include "blt_common/snp_pos_info.hh"
#include "strelka/strelka_shared.hh"

#include <iosfwd>
#include <memory>


struct strelka_pile_caller {

    strelka_pile_caller(strelka_options& opt,
                        std::ostream& os);

    void
    call(const unsigned pos,
         snp_pos_info& norm_pi,
         snp_pos_info& tumor_pi);

private:
    strelka_options& _opt;
    std::auto_ptr<strelka_deriv_options> _dopt_ptr;
    std::ostream& _os;
};

#endif

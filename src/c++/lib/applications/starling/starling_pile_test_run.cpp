//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2018 Illumina, Inc.
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
/// note coding convention for all ranges '_pos fields' is:
/// XXX_begin_pos is zero-indexed position at the beginning of the range
/// XXX_end_pos is zero-index position 1 step after the end of the range
///


#include "starling_pile_test_run.hh"

#include "gvcf_locus_info.hh"
#include "starling_common/PileupCleaner.hh"

#include <cassert>
#include <cctype>

#include <iostream>



starling_pile_caller::
starling_pile_caller(starling_options& opt,
                     std::ostream& os)
    : _opt(opt)
    , _os(os)
{
    // set default parameters:
    opt.bsnp_diploid_theta = 0.001;

    static const std::string ref_seq("ACGT");
    reference_contig_segment ref;
    ref.seq()=ref_seq;
    _dopt_ptr.reset(new starling_deriv_options(opt));
}



void
starling_pile_caller::
call(
    const unsigned /*pos*/,
    const snp_pos_info& pi)
{
    static PileupCleaner pileupCleaner(_opt);
    // recreate data cache:
    CleanedPileup cpi;

    const bool is_include_tier2(false);
    pileupCleaner.CleanPileup(pi,is_include_tier2,cpi);

    //pileupCleaner.CleanPileupErrorProb(cleanedPileup);

//    const snp_pos_info& good_pi(cleanedPileup.cleanedPileup());
    const extended_pos_info& good_epi(cpi.getExtendedPosInfo());

    diploid_genotype dgt;
    _dopt_ptr->pdcaller().position_snp_call_pprob_digt(
        _opt,good_epi, dgt, _opt.is_all_sites());

    _os << dgt;
}

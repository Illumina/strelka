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

/// \file
/// \author Chris Saunders
///

#pragma once

#include "blt_common/blt_shared.hh"
#include "blt_common/nonref_test_call.hh"
#include "blt_common/snp_pos_info.hh"

#include <iosfwd>


namespace NR2TEST
{
enum index_t
{
    REF,             // reference allele only
    NONREF_MF,       // reference allele mixed with single signal allele
    NONREF_MF_NOISE, // reference allele mixed with site specific error
    NONREF_OTHER,    // reference allele mixed with non-signal allele (lhood is approximated to 0 for this state)
    SIZE
};

inline
const char*
label(const index_t i)
{
    switch (i)
    {
    case REF:
        return "ref";
    case NONREF_MF:
        return "nonref";
    case NONREF_MF_NOISE:
        return "noise";
    case NONREF_OTHER:
        return "nonref-other";
    default:
        return "xxx";
    }
}
}



/// \brief Call a snp at a position under the assumption that up to
/// one nonref allele could occur at any frequency. When a snp is
/// found, optionally also provide the allele MLEs.
///
/// Method is optimized to pick up low-frequency snps, but will also
/// inevitably call higher frequency germline variation (without
/// correct priors) if present.
///
/// \param nonref_variant_rate The expected non-reference variant frequency (suggested start point: 0.000001)
/// \param min_nonref_freq The minimum non-reference allele frequency considered
/// \param nonref_site_error_rate The expected rate of erroneous non-reference allele sites applied to the nonref model. At error sites a nonref allele is expected in the frequency range [0,decay_freq], with a probability that linearly decays to zero at decay_freq. (suggested start point: 0.0001)
/// \param nonref_site_error_decay_freq The decay_freq used for the site-error state as described above. (suggested start point 0.01)
///
void
position_nonref_2allele_test(
    const snp_pos_info& pi,
    const double nonref_variant_rate,
    const double min_nonref_freq,
    const double nonref_site_error_rate,
    const double nonref_site_error_decay_freq,
    const bool /*is_always_test*/,
    nonref_test_call& nrc);


void
write_nonref_2allele_test(
    const snp_pos_info& pi,
    const nonref_test_call& nrc,
    std::ostream& os);

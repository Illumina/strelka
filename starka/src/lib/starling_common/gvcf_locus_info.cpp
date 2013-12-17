// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/downloads/sequencing/licenses/>.
//

/// \file
///
/// \author Chris Saunders
///


#include "starling_common/gvcf_locus_info.hh"

#include <iostream>


void
shared_modifiers::
write_filters(std::ostream& os) const {

    if (filters.none()) {
        os << "PASS";
        return;
    }

    bool is_sep(false);
    for (unsigned i(0); i<VCF_FILTERS::SIZE; ++i) {
        if (filters.test(i)) {
            if (is_sep) { os << ";"; }
            else       { is_sep=true; }
            os << VCF_FILTERS::get_label(i);
        }
    }
}



std::ostream&
operator<<(std::ostream& os,
           const shared_modifiers& shmod) {

    os << "gqx: " << shmod.gqx
       << " gq: " << shmod.gq
       << " max_gt: " << DIGT::label(shmod.max_gt);

    os << " filters: ";
    shmod.write_filters(os);

    return os;
}



std::ostream&
operator<<(std::ostream& os,
           const site_modifiers& smod) {

    os << static_cast<shared_modifiers>(smod) << '\n';

    os << "is_unknown: " << smod.is_unknown;
    os << " is_covered: " << smod.is_covered;
    os << " is_used_coverage: " << smod.is_used_covered;
    os << " is_zero_ploidy: " << smod.is_zero_ploidy;
    os << " is_block: " << smod.is_block;

    if (smod.modified_gt != MODIFIED_SITE_GT::NONE) {
        os << " modgt: " << MODIFIED_SITE_GT::get_label(smod.modified_gt);
    }

    return os;
}

std::ostream&
operator<<(std::ostream& os,
           const site_info& si) {
    os << "pos: " << (si.pos+1) << " " << si.get_gt();
    return os;
}


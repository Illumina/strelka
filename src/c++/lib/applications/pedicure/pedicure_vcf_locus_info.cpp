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


#include "pedicure_vcf_locus_info.hh"

#include <iostream>


void
pedicure_shared_modifiers::
write_filters(
    std::ostream& os) const
{
    if (filters.none())
    {
        os << "PASS";
        return;
    }

    bool is_sep(false);
    for (unsigned i(0); i<PEDICURE_VCF_FILTERS::SIZE; ++i)
    {
        if (! filters.test(i)) continue;

        if (is_sep)
        {
            os << ";";
        }
        else
        {
            is_sep=true;
        }
        os << PEDICURE_VCF_FILTERS::get_label(i);
    }
}



std::ostream&
operator<<(
    std::ostream& os,
    const pedicure_shared_modifiers& shmod)
{
    os << " filters: ";
    shmod.write_filters(os);

    return os;
}

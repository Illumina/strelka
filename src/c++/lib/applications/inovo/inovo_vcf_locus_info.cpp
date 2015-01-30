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


#include "inovo_vcf_locus_info.hh"

#include <iostream>


void
inovo_shared_modifiers::
write_filters(
    std::ostream& os) const
{
    if (filters.none())
    {
        os << "PASS";
        return;
    }

    bool is_sep(false);
    for (unsigned i(0); i<STRELKA_VCF_FILTERS::SIZE; ++i)
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
        os << STRELKA_VCF_FILTERS::get_label(i);
    }
}

void
inovo_shared_modifiers::
write_feature(
    std::ostream& os) const
{
    os << "\n #FEAT ";
    for (auto it = _featureVal.cbegin(); it != _featureVal.cend(); ++it)
        os << STRELKA_VQSR_FEATURES::get_feature_label(it->first) << "=" << it->second << "; ";
    os << "\n";
}



std::ostream&
operator<<(
    std::ostream& os,
    const inovo_shared_modifiers& shmod)
{
    os << " filters: ";
    shmod.write_filters(os);

    return os;
}

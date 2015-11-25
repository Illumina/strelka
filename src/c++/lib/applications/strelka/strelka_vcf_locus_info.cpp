// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2016 Illumina, Inc.
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


#include "strelka_vcf_locus_info.hh"

#include <iostream>


void
strelka_shared_modifiers::
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
strelka_shared_modifiers::
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
    const strelka_shared_modifiers& shmod)
{
    os << " filters: ";
    shmod.write_filters(os);

    return os;
}

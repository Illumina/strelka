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

///
/// note coding convention for all ranges '_pos fields' is:
/// XXX_begin_pos is zero-indexed position at the beginning of the range
/// XXX_end_pos is zero-index position 1 step after the end of the range
///

#pragma once


// subclasses of this object can be passed around to specify sample
// information where required, the base-class handles the starling
// case of unlabeled single sample analysis.
//
struct SampleSetSummary
{
    SampleSetSummary() {}

    virtual
    ~SampleSetSummary() {}

    virtual
    unsigned
    size() const
    {
        return 1;
    }

    virtual
    const char*
    get_label(
        const unsigned /*i*/) const
    {
        return "";
    }

    virtual
    const char*
    get_prefix(
        const unsigned /*i*/,
        const bool /*is_tier1*/) const
    {
        return "";
    }
};

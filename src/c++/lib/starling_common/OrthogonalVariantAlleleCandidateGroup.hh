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

#pragma once

#include "starling_common/indel.hh"
#include "starling_common/IndelBuffer.hh"
#include "blt_util/known_pos_range2.hh"

#include <iosfwd>
#include <vector>


/// A generalization of overlapping alleles:
///
struct OrthogonalVariantAlleleCandidateGroup
{
    typedef IndelBuffer::const_iterator AlleleIter_t;

    /// \return The merged reference range of all alleles in the group
    known_pos_range
    getReferenceRange() const;

    const IndelKey&
    key(
        const unsigned index) const
    {
        return iter(index)->first;
    }

    const IndelData&
    data(
        const unsigned index) const
    {
        assert(index < size());
        return getIndelData(iter(index));
    }

    const AlleleIter_t&
    iter(
        const unsigned index) const
    {
        assert(index < size());
        return alleles[index];
    }

    unsigned
    size() const
    {
        return alleles.size();
    }

    bool
    empty() const
    {
        return alleles.empty();
    }

    void
    clear()
    {
        alleles.clear();
    }

    void
    addVariantAllele(
        const AlleleIter_t alleleIter)
    {
        alleles.push_back(alleleIter);
    }

    std::vector<AlleleIter_t> alleles;
};

std::ostream& operator<<(
    std::ostream& os, const OrthogonalVariantAlleleCandidateGroup& group);

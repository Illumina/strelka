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

#pragma once

/// first push into a generalized Genotype which
/// can handle N overlapping alleles, plus consistent
/// representation of the refernece as 'just another allele'
struct AlleleGroupGenotype
{
    static const uint8_t MAX_GENOTYPE_COUNT = 6;

    bool
    isVariantIndel() const
    {
        return (variantAlleleQuality != 0);
    }

    unsigned maxGenotypeIndex;
    unsigned maxGenotypeIndexPolymorphic;

    double variantAlleleQuality;
    double genotypeQuality;
    double genotypeQualityPolymorphic;
    double posteriorProb[MAX_GENOTYPE_COUNT];
    double phredLoghood[MAX_GENOTYPE_COUNT];
};

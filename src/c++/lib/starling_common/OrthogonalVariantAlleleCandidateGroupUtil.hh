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

#include "OrthogonalVariantAlleleCandidateGroup.hh"


/// find all readIds for which a likelihood has been computed for at least one allele in this group
///
void
getAlleleGroupUnionReadIds(
    const unsigned sampleId,
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    std::set<unsigned>& readIds,
    const bool isTier1Only);


/// find all readIds for which a likelihood has been computed for all alleles in this group
///
void
getAlleleGroupIntersectionReadIds(
    const unsigned sampleId,
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    std::set<unsigned>& readIds,
    const bool isTier1Only);


/// enumerate read likelihood P(read | allele) for read 'readId' over all alleles in 'alleleGroup'
///
/// \param lhood[out] (normalized) likelihood for each allele, set to dimension "alleleGroup.alleles.size() + 1",
///                   with an extra reference allele state represented at the end of the array
///
void
getAlleleLikelihoodsFromRead(
    const unsigned sampleId,
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    const unsigned readId,
    std::vector<double>& lhood);


/// ranks the alleles in the input 'alleles' set according
/// to supporting read evidence in sample 'sampleId'
///
/// \param[in] sampleId index of the sample from which ranking evidence will be drawn
/// \param[in] alleleGroup unsorted list of input alleles
///
void
rankOrthogonalAllelesInSample(
    const unsigned sampleId,
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    OrthogonalVariantAlleleCandidateGroup& rankedAlleleGroup,
    unsigned& referenceRank);


/// perform the ranking as above, then select the top
/// N alleles (including the reference) and return
/// the top N (or N-1) non-reference alleles
void
selectTopOrthogonalAllelesInSample(
    const unsigned sampleId,
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    const unsigned selectionSize,
    OrthogonalVariantAlleleCandidateGroup& topAlleleGroup);


/// augment alleleGroup with overlapping alleles that have position other than 'pos',
/// then rerank and re-select top groupLocusPloidy alleles
///
/// \return true if every alt allele which otherwise qualifies was included based on
///         forming an orthogonal clique with the full allele set
///
bool
addAllelesAtOtherPositions(
    const pos_t pos,
    const pos_t largest_total_indel_ref_span_per_read,
    const unsigned sampleIndex,
    const int callerPloidy,
    const IndelBuffer& indelBuffer,
    OrthogonalVariantAlleleCandidateGroup& alleleGroup);

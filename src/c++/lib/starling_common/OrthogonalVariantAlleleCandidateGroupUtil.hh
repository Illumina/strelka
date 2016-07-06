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


/// refine alleleGroup to consider alts shared by all alleles, then rerank
/// and re-select top groupLocusPloidy alleles
///
void addAllelesAtOtherPositions(
    const pos_t pos,
    const pos_t largest_total_indel_ref_span_per_read,
    const unsigned sampleId,
    const int groupLocusPloidy,
    const IndelBuffer& indelBuffer,
    OrthogonalVariantAlleleCandidateGroup& alleleGroup);

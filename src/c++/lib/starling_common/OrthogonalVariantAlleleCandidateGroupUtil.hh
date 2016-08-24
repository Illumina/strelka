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
    const unsigned sampleIndex,
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    std::set<unsigned>& readIds,
    const bool isTier1Only);


/// find all readIds for which a likelihood has been computed for all alleles in this group
///
void
getAlleleGroupIntersectionReadIds(
    const unsigned sampleIndex,
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    std::set<unsigned>& readIds,
    const bool isTier1Only);


/// find set of read ids which support the entire set of alleles in alleleGroup
///
void
getAlleleGroupSupportingReadIds(
    const unsigned sampleIndex,
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    std::set<unsigned>& readIds,
    const bool isTier1Only);


/// enumerate (log of) read likelihood P(read | allele) for read 'readId' over all ref + all alt alleles in 'alleleGroup'
///
/// \param alleleLogLhood[out] log likelihood for each allele, set to dimension "alleleGroup.alleles.size() + 1",
///                   with an extra reference allele state represented at the end of the array
///
void
getAlleleLogLhoodFromRead(
    const unsigned sampleIndex,
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    const unsigned readId,
    std::vector<double>& alleleLogLhood);


/// same as above except allele likelihoods are is normalized
///
void
getAlleleNaivePosteriorFromRead(
    const unsigned sampleIndex,
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    const unsigned readId,
    std::vector<double>& alleleProb);


/// ranks the alleles in the input 'alleles' set according
/// to supporting read evidence in sample 'sampleId'
///
/// \param[in] sampleIndex index of the sample from which ranking evidence will be drawn
/// \param[in] alleleGroup unsorted list of input alleles
/// \param[out] referenceRank rank of the reference allele if it did exist in alleleGroup,
///                           for instances referenceRank of 0 indicates that reference is the most likely allele
///
void
rankOrthogonalAllelesInSample(
    const unsigned sampleIndex,
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    OrthogonalVariantAlleleCandidateGroup& rankedAlleleGroup,
    unsigned& referenceRank);


/// perform the ranking as above, then select the top
/// N alleles (including the reference) and return
/// the top N (or N-1) non-reference alleles
void
selectTopOrthogonalAllelesInSample(
    const unsigned sampleIndex,
    const OrthogonalVariantAlleleCandidateGroup& inputAlleleGroup,
    const unsigned selectionSize,
    OrthogonalVariantAlleleCandidateGroup& topAlleleGroup);


/// In each sample, select the top N alleles, N = ploidy. Aggregate these
/// top alleles over all samples, and use an approximate global ranking based
/// on the within-sample rankings
void
selectTopOrthogonalAllelesInAllSamples(
    const unsigned sampleCount,
    const std::vector<unsigned>& callerPloidy,
    const OrthogonalVariantAlleleCandidateGroup& inputAlleleGroup,
    OrthogonalVariantAlleleCandidateGroup& topAlleleGroup,
    std::vector<unsigned>& topVariantAlleleIndexPerSample);


/// augment alleleGroup with overlapping alleles that have position other than 'pos',
/// then rerank and re-select top groupLocusPloidy alleles
///
/// \return true if every alt allele which otherwise qualifies was included based on
///         forming an orthogonal clique with the full allele set
///
bool
addAllelesAtOtherPositions(
    const unsigned sampleCount,
    const std::vector<unsigned>& callerPloidy,
    const pos_t pos,
    const pos_t largest_total_indel_ref_span_per_read,
    const IndelBuffer& indelBuffer,
    OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    std::vector<IndelKey>& topVariantAllelePerSample);

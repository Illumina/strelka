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

#include "OrthogonalVariantAlleleCandidateGroup.hh"
#include "blt_util/RegionTracker.hh"


//#define DEBUG_INDEL_OVERLAP


/// Find all readIds for which a likelihood has been computed for at least one allele in this group
///
void
getAlleleGroupUnionReadIds(
    const unsigned sampleIndex,
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    std::set<unsigned>& readIds,
    const bool isTier1Only);


/// Find all readIds for which a likelihood has been computed for all alleles in \p alleleGroup
///
void
getAlleleGroupIntersectionReadIds(
    const unsigned sampleIndex,
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    std::set<unsigned>& readIds,
    const bool isTier1Only,
    const unsigned minDistanceFromReadEdge = 0);


/// Find set of read ids which support the entire set of alleles in alleleGroup
///
void
getAlleleGroupSupportingReadIds(
    const unsigned sampleIndex,
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    std::set<unsigned>& readIds,
    const bool isTier1Only);


/// Enumerate (log of) read likelihood P(read | allele) for read \p readId over all ref + all alt alleles in
/// \p alleleGroup
///
/// In case of any allele for which P(read | allele) has not been computed, use P(read | ref) as an approximation.
///
/// \param[out] alleleLogLhood Log likelihood for each allele, set to dimension "alleleGroup.alleles.size() + 1",
///                            with an extra reference allele state represented at the beginning of the array. Input
///                            value of array is ignored.
void
getAlleleLogLhoodFromRead(
    const unsigned sampleIndex,
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    const unsigned readId,
    std::vector<double>& alleleLogLhood);


/// Like getAlleleLogLhoodFromRead, except that allele likelihoods are normalized
///
void
getAlleleNaivePosteriorFromRead(
    const unsigned sampleIndex,
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    const unsigned readId,
    std::vector<double>& alleleProb);


/// Rank the input allele set according to supporting read evidence in the indicated sample.
///
/// \param[in] sampleIndex Zero-based index of the sample from which ranking evidence will be drawn
/// \param[in] alleleGroup Unsorted list of alternate alleles, at least one alternate allele must be provided
/// \param[out] rankedAlleleGroup Ranked list of the input alleles, this is object is cleared on input and has size equal to
///                               alleleGroup on output
/// \param[out] referenceRank Zero-based rank of the reference allele if it did exist in alleleGroup, for instance
///                           \p referenceRank of 0 indicates that reference is the most likely allele.
///
void
rankOrthogonalAllelesInSample(
    const unsigned sampleIndex,
    const OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    OrthogonalVariantAlleleCandidateGroup& rankedAlleleGroup,
    unsigned& referenceRank);


/// Rank alleles as describe in rankOrthogonalAllelesInSample, then select the top N (\p selectionSize) alleles
/// (including the reference) and return the top N (or N-1) non-reference alleles
///
/// \param[in] sampleIndex Zero-based index of the sample from which ranking evidence will be drawn
/// \param[in] inputAlleleGroup Unsorted list of input alleles, at least one alternate allele must be provided
/// \param[in] selectionSize Number of top-ranking alleles to select in this sample
/// \param[out] topAlleleGroup Ranked subset of the input alleles. This is cleared on input. On output its size will be
///                            \p selectionSize or less.
void
selectTopOrthogonalAllelesInSample(
    const unsigned sampleIndex,
    const OrthogonalVariantAlleleCandidateGroup& inputAlleleGroup,
    const unsigned selectionSize,
    OrthogonalVariantAlleleCandidateGroup& topAlleleGroup);


/// Select the top N alleles in each sample, N = ploidy. Aggregate these top alleles over all samples, and use an
/// approximate global ranking based on the within-sample rankings
///
/// \param[in] callerPloidyPerSample Array of size \p sampleCount, providing the ploidy to be used for calling each sample
/// \param[in] inputAlleleGroup Unsorted alt alleles from all samples, at least one alt allele must be provided
/// \param[out] topAlleleGroup Top ranked allele output. This will be a sorted subset of \p inputAlleleGroup
/// \param[out] topVariantAlleleIndexPerSample Index of most likely alt per sample, where vector position reflects
/// sample index and vector value references \p topAlleleGroup order index. Object oes not need to be initialized on
/// input, it will be resized to \p sampleCount in this function.
void
selectTopOrthogonalAllelesInAllSamples(
    const unsigned sampleCount,
    const std::vector<unsigned>& callerPloidyPerSample,
    const OrthogonalVariantAlleleCandidateGroup& inputAlleleGroup,
    OrthogonalVariantAlleleCandidateGroup& topAlleleGroup,
    std::vector<unsigned>& topVariantAlleleIndexPerSample);


/// Augment alleleGroup with overlapping alleles that have position other than \p pos, then re-rank and re-select top
/// groupLocusPloidy alleles
///
/// \param[in,out] alleleGroup On input, contains filtered top-ranking alleles over all samples indexed at pos. On
/// output this is potentially augmented with additional alleles from other positions
/// \param[out] topVariantAlleleIndexPerSample Index of most likely alt per sample, where vector position reflects
/// sample index and vector value references \p topAlleleGroup order index. Object oes not need to be initialized on
/// input, it will be resized to sampleCount in this function.
///
/// \return False if any alt alleles were excluded because they did not form an orthogonal clique with
///         the full allele set
bool
addAllelesAtOtherPositions(
    const reference_contig_segment& ref,
    const unsigned sampleCount,
    const std::vector<unsigned>& callerPloidy,
    const pos_t pos,
    const pos_t largestTotalIndelRefSpanPerRead,
    const IndelBuffer& indelBuffer,
    OrthogonalVariantAlleleCandidateGroup& alleleGroup,
    std::vector<unsigned>& topVariantAlleleIndexPerSample);

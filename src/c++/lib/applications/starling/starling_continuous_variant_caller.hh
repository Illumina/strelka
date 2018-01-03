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


struct starling_continuous_variant_caller
{
    /// Get the phred-scaled p-value for the hypothesis that 'allele' was generated as sequencing error under
    /// a simple Poisson error model
    ///
    /// \param[in] alleleObservationCount Observation count of the allele in question
    /// \param[in] totalObservationCount Observation count of all alleles
    /// \param[in] expectedObservationQscore Approximate that all observations have the same error probability given by
    ///                                       this value (expressed as a phred-scaled quality score)
    ///
    /// \return The above-described phred-scaled p-value
    static
    int
    getAlleleSequencingErrorQscore(
        const unsigned alleleObservationCount,
        const unsigned totalObservationCount,
        const int expectedObservationQscore,
        const int maxQScore);

    /// Return log likelihood ratio of the variants on either strand over both strands
    static
    double
    strandBias(
        unsigned fwdAlt,
        unsigned revAlt,
        unsigned fwdOther,
        unsigned revOther);
};

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

/*
 *  \author Morten Kallberg
 */

#pragma once

#include "gvcf_locus_info.hh"
#include "starling_shared.hh"
#include "calibration/VariantScoringModelServer.hh"


/// handles site and indel filter labeling OR EVS scoring and filtering
///
/// EVS model is specified in a json container expressing parameters for
/// potentially several model types.
///
/// When EVS cannot be applied, a fallback set of hard filters are applied
/// to the variant.
///
/// Also supports a legacy EVS model expressing both a logistic regression
/// model and a rule based filter, this model is deprecated.
///
struct ScoringModelManager
{
    ScoringModelManager(
        const starling_options& opt,
        const gvcf_deriv_options& gvcfDerivedOptions);

    /// the current chromosome must be specified before handling any classifications:
    void
    resetChrom(const std::string& chrom);

    void
    classify_site(
        GermlineDiploidSiteLocusInfo& locus) const;

    void
    classify_indel(
        GermlineDiploidIndelLocusInfo& locus) const;

    void
    applyDepthFilter(
        GermlineSiteLocusInfo& locus) const;

    void
    applyDepthFilter(
        GermlineIndelLocusInfo& locus) const;

    /// simple hard-cutoff filtration rules applied to site locus in one sample
    void
    default_classify_site(
        const unsigned sampleIndex,
        const unsigned allSampleLocusDepth,
        GermlineSiteLocusInfo& locus) const;

    /// simple hard-cutoff filtration rules applied to entire site locus
    void default_classify_site_locus(
        GermlineSiteLocusInfo& locus) const;

    /// simple hard-cutoff filtration rules applied to indel locus in one sample
    void
    default_classify_indel(
        const unsigned sampleIndex,
        const unsigned allSampleLocusDepth,
        GermlineIndelLocusInfo& locus) const;

    /// simple hard-cutoff filtration rules applied to entire indel locus
    void
    default_classify_indel_locus(
        GermlineIndelLocusInfo& locus) const;

    bool
    isEVSSiteModel() const
    {
        return static_cast<bool>(_snvScoringModelPtr);
    }

    bool
    isEVSIndelModel() const
    {
        return static_cast<bool>(_indelScoringModelPtr);
    }

    double
    getMaxDepth() const
    {
        assert(isChromSet());
        return _maxChromDepth;
    }

private:
    bool
    isChromSet() const
    {
        return (not _chromName.empty());
    }

    double
    snvEVSThreshold() const
    {
        return _snvScoringModelPtr->scoreFilterThreshold();
    }

    double
    indelEVSThreshold() const
    {
        return _indelScoringModelPtr->scoreFilterThreshold();
    }

    // for setting the vcf header filters
    const gvcf_options& _opt;
    const gvcf_deriv_options& _dopt;
    bool _isReportEVSFeatures;
    bool _isRNA;

    std::string _chromName;
    double _normChromDepth = 0.;
    double _maxChromDepth = 0.;

    std::unique_ptr<VariantScoringModelServer> _snvScoringModelPtr;
    std::unique_ptr<VariantScoringModelServer> _indelScoringModelPtr;
};

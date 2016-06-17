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

/*
 *  \author Morten Kallberg
 */

#pragma once

#include "gvcf_locus_info.hh"
#include "gvcf_options.hh"

#include <vector>
#include <map>

namespace LEGACY_CALIBRATION_MODEL
{

enum var_case
{
    HetSNP,
    HomSNP,
    HetAltSNP,
    HetIns,
    HomIns,
    HetAltIns,
    HetDel,
    HomDel,
    HetAltDel,
    SIZE
};

inline
const char*
get_label(const unsigned idx)
{
    switch (idx)
    {
    case HetSNP:
        return "snphet";
    case HomSNP:
        return "snphom";
    case HetAltSNP:
        return "snphetalt";
    case HetIns:
        return "inshet";
    case HomIns:
        return "inshom";
    case HetAltIns:
        return "inshetalt";
    case HetDel:
        return "delhet";
    case HomDel:
        return "delhom";
    case HetAltDel:
        return "delhetalt";
    default:
        assert(0);
        return NULL;
    }
}

inline
const char*
get_label_header(const unsigned idx)
{
    switch (idx)
    {
    case HetSNP:
        return "het SNP";
    case HomSNP:
        return "hom SNP";
    case HetAltSNP:
        return "het-alt SNP";
    case HetIns:
        return "het insertion";
    case HomIns:
        return "hom insertion";
    case HetAltIns:
        return "het-alt insertion";
    case HetDel:
        return "het deletion";
    case HomDel:
        return "hom deletion";
    case HetAltDel:
        return "het-alt deletion";
    default:
        assert(0);
        return NULL;
    }
}

inline
GERMLINE_VARIANT_VCF_FILTERS::index_t
get_Qscore_filter(const unsigned var_case)
{
    switch (var_case)
    {
    case HetSNP:
        return GERMLINE_VARIANT_VCF_FILTERS::LowQscoreHetSNP;
    case HomSNP:
        return GERMLINE_VARIANT_VCF_FILTERS::LowQscoreHomSNP;
    case HetAltSNP:
        return GERMLINE_VARIANT_VCF_FILTERS::LowQscoreHetAltSNP;
    case HetIns:
        return GERMLINE_VARIANT_VCF_FILTERS::LowQscoreHetIns;
    case HomIns:
        return GERMLINE_VARIANT_VCF_FILTERS::LowQscoreHomIns;
    case HetAltIns:
        return GERMLINE_VARIANT_VCF_FILTERS::LowQscoreHetAltIns;
    case HetDel:
        return GERMLINE_VARIANT_VCF_FILTERS::LowQscoreHetDel;
    case HomDel:
        return GERMLINE_VARIANT_VCF_FILTERS::LowQscoreHomDel;
    case HetAltDel:
        return GERMLINE_VARIANT_VCF_FILTERS::LowQscoreHetAltDel;
    default:
        assert(0);
        return GERMLINE_VARIANT_VCF_FILTERS::LowGQX;
    }
}
}


typedef std::map<std::string, double> featuremap;

/// this captures a number of concepts, including a logistic regression model
/// a rule based filtration service and a variant object to model features
/// translation
///
struct LogisticAndRuleScoringModels
{
    LogisticAndRuleScoringModels(
        const gvcf_deriv_options& dopt,
        const std::string& model_file,
        const std::string& name);

    void
    score_site_instance(
        const GermlineDiploidSiteCallInfo& si,
        GermlineDiploidSiteSimpleGenotypeInfo& smod) const;

    void
    score_indel_instance(
        const GermlineDiploidIndelCallInfo& ii,
        GermlineDiploidIndelSimpleGenotypeInfo& call) const;

    bool is_logistic_model() const;

    int
    get_var_threshold(
        const LEGACY_CALIBRATION_MODEL::var_case& my_case) const;

private:

    typedef std::map<std::string, std::map<std::string, featuremap > > parmap;

    double normal_depth() const;

    int
    logistic_score(
        const LEGACY_CALIBRATION_MODEL::var_case var_case,
        const featuremap& features) const;

    void
    do_site_rule_model(
        const featuremap& cutoffs,
        const GermlineDiploidSiteCallInfo& si,
        GermlineDiploidSiteSimpleGenotypeInfo& smod) const;

    void
    do_indel_rule_model(
        const featuremap& cutoffs,
        const GermlineDiploidIndelCallInfo& ii,
        GermlineDiploidIndelSimpleGenotypeInfo& call) const;

    /// Transform the features with the specified scaling parameters that were used to standardize
    /// the dataset to zero mean and unit variance: newVal = (oldVal-centerVal)/scaleVal.
    featuremap
    normalize(
        const featuremap& features,
        const featuremap& adjust_factor,
        const featuremap& norm_factor) const;

    double
    log_odds(
        const featuremap& features,
        const featuremap& coeffs) const;

    void
    apply_site_qscore_filters(
        const LEGACY_CALIBRATION_MODEL::var_case my_case,
        const GermlineDiploidSiteCallInfo& si,
        GermlineDiploidSiteSimpleGenotypeInfo& smod) const;

    void
    apply_indel_qscore_filters(
        const LEGACY_CALIBRATION_MODEL::var_case my_case,
        const GermlineDiploidIndelCallInfo& ii,
        GermlineDiploidIndelSimpleGenotypeInfo& call) const;

    const gvcf_deriv_options& _dopt;
    std::string _modelType;
    parmap _pars;
};

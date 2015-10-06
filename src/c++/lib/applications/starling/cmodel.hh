// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Starka
// Copyright (c) 2009-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//
/*
 *  Created on: Jan 15, 2014
 *      Author: Morten Kallberg
 */

#pragma once

#include "gvcf_locus_info.hh"
#include "gvcf_options.hh"

#include <vector>
#include <map>

namespace CALIBRATION_MODEL
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
VCF_FILTERS::index_t
get_Qscore_filter(const unsigned var_case)
{
    switch (var_case)
    {
    case HetSNP:
        return VCF_FILTERS::LowQscoreHetSNP;
    case HomSNP:
        return VCF_FILTERS::LowQscoreHomSNP;
    case HetAltSNP:
        return VCF_FILTERS::LowQscoreHetAltSNP;
    case HetIns:
        return VCF_FILTERS::LowQscoreHetIns;
    case HomIns:
        return VCF_FILTERS::LowQscoreHomIns;
    case HetAltIns:
        return VCF_FILTERS::LowQscoreHetAltIns;
    case HetDel:
        return VCF_FILTERS::LowQscoreHetDel;
    case HomDel:
        return VCF_FILTERS::LowQscoreHomDel;
    case HetAltDel:
        return VCF_FILTERS::LowQscoreHetAltDel;
    default:
        assert(0);
        return VCF_FILTERS::LowGQX;
    }
}
}

typedef std::map<std::string, double> featuremap;
typedef std::map<std::string, std::map<std::string, featuremap > > parmap;

struct c_model
{
    c_model(
        const std::string& name,
        const std::string& type,
        const gvcf_deriv_options& init_dopt) :
        dopt(init_dopt),
        model_name(name),
        model_type(type)
    {}

    void
    add_parameters(const parmap& myPars);

    void
    score_site_instance(
        const digt_site_info& si,
        digt_call_info& smod) const;

    void
    score_indel_instance(
        const digt_indel_info& ii,
        digt_indel_call& call) const;

    bool is_logistic_model() const;

    int
    get_var_threshold(
        const CALIBRATION_MODEL::var_case& my_case) const;

private:

    double normal_depth() const;

    int
    logistic_score(
        const CALIBRATION_MODEL::var_case var_case,
        const featuremap& features) const;

    void
    do_site_rule_model(
        const featuremap& cutoffs,
        const digt_site_info& si,
        digt_call_info& smod) const;

    void
    do_indel_rule_model(
        const featuremap& cutoffs,
        const digt_indel_info& ii,
        digt_indel_call& call) const;

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
        const CALIBRATION_MODEL::var_case my_case,
        const digt_site_info& si,
        digt_call_info& smod) const;

    void
    apply_indel_qscore_filters(
        const CALIBRATION_MODEL::var_case my_case,
        const digt_indel_info& ii,
        digt_indel_call& call) const;

    const gvcf_deriv_options& dopt;
    std::string model_name;
    std::string model_type;
    parmap pars;
};

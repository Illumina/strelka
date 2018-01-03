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

/// \file
/// \author Chris Saunders
///

#pragma once

#include "gvcf_options.hh"
#include "starling_common/starling_base_shared.hh"


struct starling_options : public starling_base_options
{
    starling_options()
    {
        bsnp_ssd_no_mismatch = 0.35;
        bsnp_ssd_one_mismatch = 0.6;
        mismatchDensityFilterMaxMismatchCount = 2;
        mismatchDensityFilterFlankSize = 20;
        is_min_vexp = true;
        min_vexp = 0.25;

        // In practice, we find that increasing the indel -> ref error rate relative to baseline improves performance
        // of the germline model. This is likely due to various forms of reference bias in the realignment and
        // scoring processes.
        //
        // Setting this factor above 1 acts to compensate against this reference allele bias.
        //
        isIndelRefErrorFactor = true;
        indelRefErrorFactor = 1.8;

        // turn on short haplotyping
        isHaplotypingEnabled = true;

        // the germline indel error model defaults to a static version of the adaptive indel error estimation values
        indel_error_model_name = "adaptiveDefault";
    }

    bool
    is_bsnp_diploid() const override
    {
        return is_ploidy_prior;
    }

    bool
    is_all_sites() const override
    {
        return true;
    }

    bool
    is_compute_germline_scoring_metrics() const override
    {
        return (isReportEVSFeatures || (! snv_scoring_model_filename.empty()) || (! indel_scoring_model_filename.empty()));
    }

    const AlignmentFileOptions&
    getAlignmentFileOptions() const override
    {
        return alignFileOpt;
    }

    AlignmentFileOptions alignFileOpt;

    // empirical scoring models
    std::string snv_scoring_model_filename;
    std::string indel_scoring_model_filename;

    /// \brief Turn on read backed variant phasing if true
    bool isUseVariantPhaser = false;

    /// \brief Apply special behaviors for RNA-Seq analysis if true
    bool isRNA = false;

    /// This is the absolute value limit of the likelihood-based strand bias score range. It must be a positive value.
    const double maxAbsSampleVariantStrandBias = 99;

    gvcf_options gvcf;
};


/// data deterministically derived from the input options:
///
struct starling_deriv_options : public starling_base_deriv_options
{
    typedef starling_base_deriv_options base_t;

    explicit
    starling_deriv_options(const starling_options& opt)
        : base_t(opt),
          gvcf(opt.gvcf, opt.isRNA)
    {}

    gvcf_deriv_options gvcf;
};

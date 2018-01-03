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

#include "blt_util/blt_types.hh"
#include "blt_util/PolymorphicObject.hh"
#include "blt_util/pos_range.hh"
#include "blt_util/seq_util.hh"

#include <cassert>

#include <memory>
#include <vector>

extern const unsigned MAX_FLANK_SIZE;


namespace LOG_LEVEL
{
enum index_t
{
    DEFAULT = 0, ///< errors and low-frequency warnings
    ALLWARN = 1  ///< all other warnings
};
}



struct blt_options : public PolymorphicObject
{
    void
    validate() const
    {
        // this should never be true once the object is fully constructed, TODO: can we move this into a state enum so it *can't* be true?:
        assert(! (is_compute_germline_scoring_metrics() && is_compute_somatic_scoring_metrics));

        assert((hetVariantFrequencyExtension >= 0.) && (hetVariantFrequencyExtension < 0.5));
    }

    virtual
    bool
    is_compute_germline_scoring_metrics() const
    {
        return false;
    }

    virtual
    bool
    is_bsnp_diploid() const
    {
        return false;
    }

    bool
    is_dependent_eprob() const
    {
        return (is_bsnp_diploid() &&
                (bsnp_ssd_no_mismatch>0. || bsnp_ssd_one_mismatch>0));
    }

    double bsnp_diploid_theta = 0.001;
    double bsnp_ssd_no_mismatch = 0;
    double bsnp_ssd_one_mismatch = 0;

    /// Expand the variant frequency range over which heterozygous variants are modeled. When this is 0 the standard
    /// DNA-seq caller behavior of modeling hets as a point process at 0.5 is used. This can be useful for representing
    /// RNA-seq ASE and other factors.
    ///
    /// This argument must be in range [0,0.5)
    double hetVariantFrequencyExtension = 0;

    bool
    isHetVariantFrequencyExtensionDefined() const
    {
        return (hetVariantFrequencyExtension>0);
    }

    /// When hetVariantFrequencyExtension is defined, the continuous frequency range for heterozygous variants is
    /// represented by a frequency grid with a grid increment distance less than this value
    ///
    /// The additional "fudge factor" added to 0.05 is used because the requested hetVariantFrequencyExtension is
    /// often a multiple of 0.05, which would force us the use the finest grid spacing (0.025) in nearly all
    /// normal cases otherwise.
    const double maxHetVariantFrequencyIncrement = 0.0501;

    int minBasecallErrorPhredProb = 17;
    int minMappingErrorPhredProb = 20;

    /// Don't use a basecall during SNV calling if mismatch count>max_win_mismatch
    /// within a window of max_win_mismatch_flank_size flanking bases
    unsigned mismatchDensityFilterMaxMismatchCount = 0;
    unsigned mismatchDensityFilterFlankSize = 0;

    bool
    isMismatchDensityFilter() const
    {
        return (mismatchDensityFilterFlankSize > 0);
    }

    /// If true, use reads with unmapped mates for variant calling
    bool includeSingletonReads = false;

    /// If true, use non proper-pair reads for variant calling
    bool includeAnomalousReads = false;

    bool is_min_vexp = false;
    double min_vexp = 0;

    /// Setting to change warning sensitivity
    LOG_LEVEL::index_t verbosity = LOG_LEVEL::DEFAULT;

    std::string cmdline;

    /// If true, then reads above max_input_depth are filtered out (option used mostly for amplicon calling)
    bool is_max_input_depth = false;

    /// If is_max_input_depth, then filter out reads above this depth
    unsigned max_input_depth = 0;

    /// If true, then report production and development EVS features in VCF output files
    bool isReportEVSFeatures = false;
    bool is_compute_somatic_scoring_metrics = false;
};



struct pprob_digt_caller;


/// data deterministically derived from the user input options
struct blt_deriv_options
{
    blt_deriv_options(
        const blt_options& opt);

    ~blt_deriv_options();

    const pprob_digt_caller&
    pdcaller() const
    {
        return *(_pdcaller.get());
    }

private:
    std::unique_ptr<pprob_digt_caller> _pdcaller;
};



struct blt_read_counts
{
    void
    report(std::ostream& os) const;

    unsigned subsample_filter = 0;
    unsigned primary_filter = 0;
    unsigned duplicate = 0;
    unsigned unmapped = 0;
    unsigned secondary = 0;
    unsigned supplement = 0;
    unsigned large_ref_deletion = 0;
    unsigned align_score_filter = 0;

    // Floating means the read is indicated as mapped, but has no "M"
    // in the cigar string. Typically inside of an insertion.
    unsigned floating = 0;

    // if optional setting is given to filter out reads once a certain depth
    // is exceeded, the number of reads filtered are enumerated here:
    unsigned max_depth = 0;
    unsigned used = 0;
};

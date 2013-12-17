// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

/// \file
///
/// \author Chris Saunders
///

#include "blt_common/blt_arg_validate.hh"
#include "blt_util/log.hh"

#include <cstdlib>

#include <iostream>
#include <iomanip>
#include <sstream>


const double MAX_MONOPLOID_THETA(1.);
const double MAX_DIPLOID_THETA(0.38068); // solution of: 0 = 1/3*(theta**2) + (5/2)*theta - 1



void
check_option_arg_range(const prog_info& pinfo,
                       const double val,
                       const char* label,
                       const double min,
                       const double max) {

    if ((val >= min) && (val <= max)) return;

    std::ostringstream oss;
    oss << std::setprecision(10);
    oss << "Value provided for '" << label << "': '" << val
        << "', is not in expected range: [ " << min << " , " << max << " ]";
    pinfo.usage(oss.str().c_str());
}



void
validate_blt_opt(const prog_info& pinfo,
                 const blt_options& opt) {

    if (opt.is_min_win_qscore && opt.min_win_qscore_flank_size > MAX_FLANK_SIZE) {
        std::ostringstream oss;
        oss << "min-window-qscore flank size exceeds max value of: " << MAX_FLANK_SIZE;
        pinfo.usage(oss.str().c_str());
    }

    if (opt.is_max_win_mismatch && opt.max_win_mismatch_flank_size > MAX_FLANK_SIZE) {
        std::ostringstream oss;
        oss << "max-window-mismatch flank size exceeds max value of: " << MAX_FLANK_SIZE;
        pinfo.usage(oss.str().c_str());
    }

    if (opt.is_adis_win_lrt && opt.adis_win_lrt_flank_size > MAX_FLANK_SIZE) {
        std::ostringstream oss;
        oss << "anom-distro-window-lrt flank size exceeds max value of: " << MAX_FLANK_SIZE;
        pinfo.usage(oss.str().c_str());
    }

    if (opt.bsnp_monoploid_theta>MAX_MONOPLOID_THETA) {
        std::ostringstream oss;
        oss << "monoploid heterozygosity exceeds maximum value of: " << MAX_MONOPLOID_THETA;
        pinfo.usage(oss.str().c_str());
    }

    if (opt.bsnp_diploid_theta>MAX_DIPLOID_THETA) {
        std::ostringstream oss;
        oss << "diploid heterozygosity exceeds maximum value of: " << MAX_DIPLOID_THETA;
        pinfo.usage(oss.str().c_str());
    }

    if (opt.is_bsnp_nploid && opt.bsnp_nploid_ploidy <= 0) {
        pinfo.usage("ERROR:: ploidy argument to -bsnp-nploid must be greater than 0\n");
    }

    if (opt.is_bsnp_diploid_het_bias) {
        if (! opt.is_bsnp_diploid()) {
            pinfo.usage("-bsnp-diploid-het-bias does not make sense when bsnp-diploid model is not in use\n");
        }
        if ((opt.bsnp_diploid_het_bias < 0.) || (opt.bsnp_diploid_het_bias >= 0.5)) {
            pinfo.usage("-bsnp-diploid-het-bias argument must be in range [0,0.5)\n");
        }
    }

    if (opt.is_nonref_test()) {
        check_option_arg_range(pinfo,opt.nonref_variant_rate,"nonref-variant-rate",0.,1.);
        check_option_arg_range(pinfo,opt.min_nonref_freq,"min-nonref-freq",0.,1.);
    }

    const pos_range& rr(opt.user_report_range);

    if (opt.is_report_range_ref) {
        if       (! opt.is_ref_set()) {
            pinfo.usage("a reference must be specified when using -report-range-reference flag");
        } else if (rr.is_begin_pos) {
            pinfo.usage("-report-range-begin cannot be combined with -report-range-reference flag");
        } else if (rr.is_end_pos) {
            pinfo.usage("-report-range-end cannot be combined with -report-range-reference flag");
        }
    }

    if (rr.is_begin_pos) {
        if (rr.begin_pos<=0) pinfo.usage("-report-range-begin argument must be > 0");
    }

    if (rr.is_end_pos) {
        if (rr.end_pos<=0) pinfo.usage("-report-range-end argument must be > 0");
        if (rr.is_begin_pos && rr.end_pos<rr.begin_pos) {
            pinfo.usage("-report-range-end argument must be >= to -report-range-begin argument");
        }
    }
}

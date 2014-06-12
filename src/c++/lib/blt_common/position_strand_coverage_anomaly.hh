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

/// \file

/// \author Chris Saunders
///
#ifndef __POSITION_STRAND_COVERAGE_ANOMALY_HH
#define __POSITION_STRAND_COVERAGE_ANOMALY_HH

#include "blt_common/snp_pos_info.hh"


/// \brief call a strand coverage anomaly with FPR = alpha
///
bool
position_strand_coverage_anomaly(const double alpha,
                                 const snp_pos_info& pi);

#if 0
/// \brief return the p-value of the strand coverage anomaly test
///
double
position_strand_coverage_anomaly_pval(const snp_pos_info& pi);
#endif

#endif

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
#ifndef __POSITION_STRAND_DISTRO_ANOMALY_HH
#define __POSITION_STRAND_DISTRO_ANOMALY_HH

#include "blt_common/snp_pos_info.hh"


double
position_strand_distro_anomaly_pval(const snp_pos_info& pi,
                                    double* ws);


/// \brief call a strand distribution anomaly test with FPR = alpha
///
/// Uses contingency table analysis with a max error prob filter
///
bool
position_strand_distro_anomaly(const double alpha,
                               const snp_pos_info& pi,
                               double* ws);

#endif

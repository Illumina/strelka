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
#ifndef __POSITION_STRAND_DISTRO_ANOMALY_LRT_HH
#define __POSITION_STRAND_DISTRO_ANOMALY_LRT_HH

#include "blt_common/blt_shared.hh"

void
position_strand_distro_anomaly_lrt_expert(const snp_pos_info& pi,
                                          double& null_loghood,
                                          double& alt_loghood,
                                          unsigned& df);

/// \brief call a strand distribution anomaly with FPR = alpha
///
bool
position_strand_distro_anomaly_lrt(const double alpha,
                                   const snp_pos_info& pi);

#endif

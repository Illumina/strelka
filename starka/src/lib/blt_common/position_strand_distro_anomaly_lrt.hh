// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2012 Illumina, Inc.
//
// This software is covered by the "Illumina Genome Analyzer Software
// License Agreement" and the "Illumina Source Code License Agreement",
// and certain third party copyright/licenses, and any user of this
// source file is bound by the terms therein (see accompanying files
// Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
// Illumina_Source_Code_License_Agreement.pdf and third party
// copyright/license notices).
//
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

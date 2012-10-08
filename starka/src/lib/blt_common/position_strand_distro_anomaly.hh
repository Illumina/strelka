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

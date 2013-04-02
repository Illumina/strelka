// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/downloads/sequencing/licenses/>.
//

/// \file

/// \author Chris Saunders
///
#ifndef __STAT_UTIL_HH
#define __STAT_UTIL_HH

bool
is_chi_sqr_reject(const double xsq,
                  const unsigned df,
                  const double alpha);

bool
is_lrt_reject_null(const double null_loghood,
                   const double alt_loghood,
                   const unsigned df,
                   const double alpha);

#endif

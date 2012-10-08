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
#include "blt_util/stat_util.hh"

#include <boost/math/distributions/chi_squared.hpp>



bool
is_chi_sqr_reject(const double xsq,
                  const unsigned df,
                  const double alpha){

    assert(xsq>=0);
    assert(df>0);

    boost::math::chi_squared dist(df);
    return ((1.-boost::math::cdf(dist,xsq)) < alpha);

#if 0
    // alternate implementation (is one faster?):
    const double xsq_crit_val(boost::math::quantile(dist,1.-alpha));
    return xsq>xsq_crit_val;
#endif
}



bool
is_lrt_reject_null(const double null_loghood,
                   const double alt_loghood,
                   const unsigned df,
                   const double alpha){

    if(df == 0) return false;
    if(null_loghood>alt_loghood) return false;

    const double log_lrt(-2.*(null_loghood-alt_loghood));

    return is_chi_sqr_reject(log_lrt,df,alpha);
}

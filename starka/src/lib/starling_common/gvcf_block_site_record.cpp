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
///
/// \author Chris Saunders
///

#include "gvcf_block_site_record.hh"
#include "blt_util/compat_util.hh"



static
bool
check_block_single_tolerance(const stream_stat& ss,
                             const int min,
                             const int tol) {
    return ((min + tol) >= ss.max());
}



static
bool
check_block_tolerance(const stream_stat& ss,
                      const double frac_tol,
                      const int abs_tol) {

    const int min(static_cast<int>(compat_round(ss.min())));
    if(check_block_single_tolerance(ss,min,abs_tol)) return true;
    const int ftol(static_cast<int>(std::floor(min * frac_tol)));
    if (ftol <= abs_tol) return false;
    return check_block_single_tolerance(ss, min, ftol);
}



static
bool
is_new_value_blockable(const bool is_new_val,
                       const int new_val,
                       const bool is_old_val,
                       const stream_stat& ss,
                       const double frac_tol,
                       const int abs_tol) {

    if(!(is_new_val && is_old_val)) return (is_new_val == is_old_val);

    stream_stat ss2(ss);
    ss2.add(new_val);
    return check_block_tolerance(ss2,frac_tol,abs_tol);
}



bool
gvcf_block_site_record::
test(const site_info& si) const {

    if(count==0) return true;
    
    // pos must be +1 from end of record:
    if((record.pos+count) != si.pos) return false;
    
    // filters must match:
    if(record.smod.filters != si.smod.filters) return false;
    
    if(0!=strcmp(record.get_gt(),si.get_gt())) return false;
    
    // coverage states must match:
    if(record.smod.is_covered != si.smod.is_covered) return false;
    if(record.smod.is_used_covered != si.smod.is_used_covered) return false;
    
    // test blocking values:
    if(! is_new_value_blockable(si.smod.is_gqx(),si.smod.gqx,
                                record.smod.is_gqx(),
                                block_gqx,frac_tol,abs_tol)) {
        return false;
    }
       
    return true;
}



void
gvcf_block_site_record::
join(const site_info& si) {
    if(count == 0) {
        record = si;
        record.smod.is_block=true;
    }

    if(si.smod.is_gqx()) {
        block_gqx.add(si.smod.gqx);
    }

    count += 1;
}



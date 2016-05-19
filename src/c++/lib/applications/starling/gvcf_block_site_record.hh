// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2016 Illumina, Inc.
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

///
/// \author Chris Saunders
///

#pragma once


#include "gvcf_locus_info.hh"
#include "gvcf_options.hh"
#include "blt_util/stream_stat.hh"


/// manages compressed site record blocks output in the gVCF
///
struct gvcf_block_site_record : public GermlineVariantSimpleGenotypeInfo
{
    explicit
    gvcf_block_site_record(const gvcf_options& opt)
        : frac_tol(static_cast<double>(opt.block_percent_tol)/100.)
        , abs_tol(opt.block_abs_tol)
    {
        reset();
    }

    void
    reset()
    {
        count=0;
        block_gqx.reset();
        block_dpu.reset();
        block_dpf.reset();
        pos=-1;
        ref=(char)0;
        gt = ".";
        has_call = is_covered = is_used_covered = is_nonref=false;
        ploidy = 0;

        GermlineVariantSimpleGenotypeInfo::clear();
    }

    /// determine if the given site could be joined to this block:
    bool test(const GermlineDiploidSiteCallInfo& si) const;

    /// add site to the current block
    void join(const GermlineDiploidSiteCallInfo& si);

    /// determine if the given site could be joined to this block:
    bool test(const GermlineContinuousSiteCallInfo& si) const;

    /// add site to the current block
    void join(const GermlineContinuousSiteCallInfo& si);

    const char* get_gt() const
    {
        return gt.c_str();
    }


    const double frac_tol;
    const int abs_tol;
    int count;
    pos_t pos;
    char ref;
    bool is_nonref;
    bool is_covered;
    bool is_used_covered;
    int ploidy;
    std::string gt;
    stream_stat block_gqx;
    stream_stat block_dpu;
    stream_stat block_dpf;

    bool has_call;
    //stream_stat _blockMQ;
};

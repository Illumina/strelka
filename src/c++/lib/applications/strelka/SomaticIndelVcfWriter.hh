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

///
/// \author Chris Saunders
///

#pragma once

#include "somatic_indel_call.hh"
#include "strelka_shared.hh"

#include "starling_common/starling_indel_report_info.hh"
#include "starling_common/starling_pos_processor_win_avg_set.hh"

#include <array>
#include <map>


/// store vcf indel record info to write out later:
struct SomaticIndelVcfInfo
{
    somatic_indel_call sindel;
    starling_indel_report_info iri;
    std::array<starling_indel_sample_report_info,2> nisri;
    std::array<starling_indel_sample_report_info,2> tisri;
};


/// delay writing indel vcf record until window data is computed:
struct SomaticIndelVcfWriter
{
    SomaticIndelVcfWriter(
        const strelka_options& opt,
        const strelka_deriv_options& dopt,
        std::ostream* osptr) :
        _opt(opt),
        _dopt(dopt),
        _osptr(osptr)
    {}

    /// return true if indel information is cached for this position
    bool
    testPos(
        const pos_t pos) const
    {
        return (_data.count(pos) != 0);
    }

    /// return true if no indel information is cache
    bool
    empty() const
    {
        return (_data.empty());
    }

    /// store an indel call
    void
    cacheIndel(
        const pos_t pos,
        const SomaticIndelVcfInfo& siInfo);

    /// add final information required
    void
    addIndelWindowData(
        const pos_t pos,
        const win_avg_set& wasNormal,
        const win_avg_set& wasTumor);

private:
    const strelka_options& _opt;
    const strelka_deriv_options& _dopt;
    std::ostream* _osptr;
    std::map<pos_t,std::vector<SomaticIndelVcfInfo>> _data;
};

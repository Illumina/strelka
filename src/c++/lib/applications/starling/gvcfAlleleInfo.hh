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

#include "germlineVariantEmpiricalScoringFeatures.hh"
#include "blt_common/position_snp_call_pprob_digt.hh"
#include "blt_util/align_path.hh"
#include "blt_util/math_util.hh"
#include "blt_util/PolymorphicObject.hh"
#include "starling_common/starling_indel_call_pprob_digt.hh"
#include "gvcf_options.hh"

#include <iosfwd>
#include <map>


/// this object contains information shared by any germline variant allele
///
/// variant here means SNV or indel
///
/// SimpleGenotype here means information pertaining to the core genotyping algorithm for SNVs or Indels (where for indels this genotype is
///    restricted to a single alt allele). This information is suplemnted with additional details requred for full VCF records in the "Call"
///    objects further below
///
/// model types where this is used include diploid/haploid and continuous single sample calling, but not contrastive (ie. tumor-normal) models
///
struct GermlineVariantAlleleInfo : public PolymorphicObject
{
    GermlineVariantAlleleInfo()
    {
        clear();
    }

    void
    clear()
    {
        gqx = 0;
        gq = 0;
        strand_bias = 0;
    }

    int gqx=0;
    int gq=0;
    double strand_bias = 0;
};


std::ostream& operator<<(std::ostream& os,const GermlineVariantAlleleInfo& shmod);


/// restrict to the case where variant is indel
struct GermlineIndelAlleleInfo : public GermlineVariantAlleleInfo
{
    GermlineIndelAlleleInfo(
        const IndelKey& indelKey,
        const IndelData& indelData,
        const starling_indel_report_info indelReportInfo,
        const starling_indel_sample_report_info& indelSampleReportInfo)
        : _indelKey(indelKey)
        , _indelData(indelData)
        , _indelReportInfo(indelReportInfo)
        , _indelSampleReportInfo(indelSampleReportInfo)
    {}

    void
    clear()
    {
        GermlineVariantAlleleInfo::clear();
        cigar.clear();
    }

    void set_hap_cigar(
        const unsigned lead=1,
        const unsigned trail=0);

    const IndelKey _indelKey;
    const IndelData _indelData;
    // TODO: make the indel overlapping code create a new call, then revert this to const
    starling_indel_report_info _indelReportInfo;
    const starling_indel_sample_report_info _indelSampleReportInfo;

    ALIGNPATH::path_t cigar;
};

std::ostream& operator<<(std::ostream& os,const GermlineIndelAlleleInfo& shi);


/// restrict to diploid calling models
struct GermlineDiploidIndelAlleleInfo : public GermlineIndelAlleleInfo
{
    GermlineDiploidIndelAlleleInfo(
        const IndelKey& indelKey,
        const IndelData& indelData,
        const starling_indel_report_info& indelReportInfo,
        const starling_indel_sample_report_info& indelSampleReportInfo,
        const GermlineDiploidIndelSimpleGenotypeInfoCore& dindel)
        : GermlineIndelAlleleInfo(indelKey, indelData, indelReportInfo, indelSampleReportInfo)
        , _dindel(dindel)
    {}

    /// TODO: this object is only here to trigger an auto copy ctor, seems like we can simplify now?
    /// TODO: Make indel_overlapper create new call objects, then revert this to const
    GermlineDiploidIndelSimpleGenotypeInfoCore _dindel;

    /// TODO: max_gt is here and in _dindel. Ugggghhh. Document the difference or get this down to one copy
    unsigned max_gt=0;
};

std::ostream& operator<<(std::ostream& os,const GermlineDiploidIndelAlleleInfo& dic);



namespace MODIFIED_SITE_GT
{

enum index_t
{
    NONE,
    UNKNOWN,
    ZERO,
    ONE
};

inline
const char*
get_label(const unsigned idx)
{
    switch (static_cast<index_t>(idx))
    {
    case ZERO:
        return "0";
    case ONE:
        return "1";
    case UNKNOWN:
        return ".";
    default:
        assert(false && "Unknown site GT value");
        return nullptr;
    }
}
}


/// restrict to the case where variant is site/SNV and calling model is diploid
struct GermlineDiploidSiteAlleleInfo : public GermlineVariantAlleleInfo
{
    GermlineDiploidSiteAlleleInfo()
    {
        clear();
    }

    void
    clear()
    {
        GermlineVariantAlleleInfo::clear();
        is_unknown=true;
        is_covered=false;
        is_used_covered=false;
        is_zero_ploidy=false;
        is_phased_region=false;
        is_phasing_insufficient_depth=false;
        modified_gt=MODIFIED_SITE_GT::NONE;
        max_gt=0;
    }

    bool
    is_gqx() const
    {
        return ((!is_unknown) && is_used_covered && (!is_zero_ploidy));
    }

    bool is_unknown;
    bool is_covered;
    bool is_used_covered;
    bool is_zero_ploidy;
    bool is_phased_region;
    bool is_phasing_insufficient_depth;

    MODIFIED_SITE_GT::index_t modified_gt;
    unsigned max_gt;
};


/// restrict to the case where variant is site/SNV and calling model is continuous
struct GermlineContinuousSiteAlleleInfo : public GermlineVariantAlleleInfo
{
    GermlineContinuousSiteAlleleInfo(unsigned totalDepth, unsigned alleleDepth, BASE_ID::index_t base)
        : _totalDepth(totalDepth)
        , _alleleDepth(alleleDepth)
        , _base(base)
    {
    }

    double variant_frequency() const
    {
        return safeFrac(_alleleDepth, _totalDepth);
    }

    unsigned _totalDepth;
    unsigned _alleleDepth;
    BASE_ID::index_t _base;
};


/// restrict to the case where variant is indel and calling model is continuous
struct GermlineContinuousIndelAlleleInfo : public GermlineIndelAlleleInfo
{
    GermlineContinuousIndelAlleleInfo(
        unsigned totalDepth,
        unsigned alleleDepth,
        const IndelKey& indelKey,
        const IndelData& indelData,
        const starling_indel_report_info& indelReportInfo,
        const starling_indel_sample_report_info& indelSampleReportInfo)
        : GermlineIndelAlleleInfo(indelKey, indelData, indelReportInfo, indelSampleReportInfo)
        , _totalDepth(totalDepth)
        , _alleleDepth(alleleDepth)
    {
        set_hap_cigar(0,0);
    }

    double variant_frequency() const
    {
        return safeFrac(_alleleDepth,_totalDepth);
    }

    unsigned _totalDepth;
    unsigned _alleleDepth;
};

std::ostream& operator<<(std::ostream& os,const GermlineDiploidSiteAlleleInfo& smod);

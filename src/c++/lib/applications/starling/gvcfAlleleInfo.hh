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


/// allele details in the context of a specific overlapping allele set
///
struct OrthogonalAlleleSetLocusReportInfoAlleleDetails
{
    std::string vcfAltSeq;
    ALIGNPATH::path_t vcfCigar;
};


/// locus summary information for an overlapping allele set which is shared between all samples
///
/// this contains information required to describe ALT and CIGAR, but not sample-specific
/// information like GT
///
struct OrthogonalAlleleSetLocusReportInfo
{
    unsigned
    getAltAlleleCount() const
    {
        return altAlleles.size();
    }

    /// 1-indexed position used to report the whole allele group
    pos_t vcfPos = 0;
    std::string vcfRefSeq;
    std::vector<OrthogonalAlleleSetLocusReportInfoAlleleDetails> altAlleles;
};



/// utility to transform indel key plus any leading or trailing match sequence
/// into an cigar alignment
void
setIndelAlleleCigar(
    const unsigned lead,
    const unsigned trail,
    const IndelKey& indelKey,
    ALIGNPATH::path_t& cigar);




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
        strand_bias = 0;
        isForcedOutput = false;
    }

    double strand_bias = 0;

    ///< this allele must appear in the VCF output:
    bool isForcedOutput = false;
};

std::ostream& operator<<(std::ostream& os,const GermlineVariantAlleleInfo& allele);


/// restrict to the case where variant is indel
struct GermlineIndelAlleleInfo : public GermlineVariantAlleleInfo
{
    GermlineIndelAlleleInfo(
        const IndelKey& initIndelKey,
        const IndelData& indelData)
        : indelKey(initIndelKey),
          indelReportInfo(indelData.getReportInfo())
    {
        isForcedOutput = indelData.isForcedOutput;
        if (indelKey.is_breakpoint())
        {
            breakpointInsertSeq = indelData.getBreakpointInsertSeq();
        }
    }

    const IndelKey indelKey;
    std::string breakpointInsertSeq;
    const AlleleReportInfo indelReportInfo;
};

std::ostream& operator<<(std::ostream& os,const GermlineIndelAlleleInfo& allele);


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


/// restrict to the case where variant is site/SNV
struct GermlineSiteAlleleInfo : public GermlineVariantAlleleInfo
{
    explicit
    GermlineSiteAlleleInfo(BASE_ID::index_t initBase)
        : base(initBase)
    {}

    BASE_ID::index_t base;
};


/// restrict to the case where variant is site/SNV and calling model is diploid
struct GermlineDiploidSiteAlleleInfo : public GermlineVariantAlleleInfo
{
    typedef GermlineVariantAlleleInfo base_t;

    GermlineDiploidSiteAlleleInfo()
    {
        clear();
    }

    void
    clear()
    {
        base_t::clear();
        is_covered=false;
        is_used_covered=false;
        is_zero_ploidy=false;
        is_phased_region=false;
        is_phasing_insufficient_depth=false;
        modified_gt=MODIFIED_SITE_GT::NONE;
        max_gt=0;
    }

    bool
    is_gqx_tmp() const
    {
        return (is_used_covered && (!is_zero_ploidy));
    }

    bool is_covered;
    bool is_used_covered;
    bool is_zero_ploidy;
    bool is_phased_region;
    bool is_phasing_insufficient_depth;

    MODIFIED_SITE_GT::index_t modified_gt;
    unsigned max_gt;
};

std::ostream& operator<<(std::ostream& os,const GermlineDiploidSiteAlleleInfo& allele);


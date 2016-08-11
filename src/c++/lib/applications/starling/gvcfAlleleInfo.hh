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
        gqx = 0;
        strand_bias = 0;
        isForcedOutput = false;
    }

    int gqx=0;
    double strand_bias = 0;

    ///< this allele must appear in the VCF output:
    bool isForcedOutput = false;
};


std::ostream& operator<<(std::ostream& os,const GermlineVariantAlleleInfo& shmod);


/// restrict to the case where variant is indel
struct GermlineIndelAlleleInfo : public GermlineVariantAlleleInfo
{
    /// new interface
    GermlineIndelAlleleInfo(
        const IndelKey& initIndelKey,
        const IndelData& indelData)
        : _indelKey(initIndelKey)
    {
        isForcedOutput = indelData.isForcedOutput;
        if (_indelKey.is_breakpoint())
        {
            breakpointInsertSeq = indelData.getBreakpointInsertSeq();
        }
    }

    /// deprecated interface
    GermlineIndelAlleleInfo(
        const IndelKey& indelKey,
        const IndelData& indelData,
        const AlleleReportInfo indelReportInfo)
        : _indelKey(indelKey)
        , _indelReportInfo(indelReportInfo)
    {
        isForcedOutput = indelData.isForcedOutput;
    }

    void
    clear()
    {
        GermlineVariantAlleleInfo::clear();
        cigar.clear();
    }

    void set_hap_cigar(
        const unsigned lead=1,
        const unsigned trail=0)
    {
        setIndelAlleleCigar(lead, trail, _indelKey, cigar);
    }

    const IndelKey _indelKey;
    std::string breakpointInsertSeq;
    // TODO: make the indel overlapping code create a new call, then revert this to const
    AlleleReportInfo _indelReportInfo;

    ALIGNPATH::path_t cigar;
};

std::ostream& operator<<(std::ostream& os,const GermlineIndelAlleleInfo& shi);


/// restrict to diploid calling models
struct GermlineDiploidIndelAlleleInfo : public GermlineIndelAlleleInfo
{
    GermlineDiploidIndelAlleleInfo(
        const IndelKey& indelKey,
        const IndelData& indelData,
        const AlleleReportInfo& indelReportInfo,
        const GermlineDiploidIndelSimpleGenotypeInfoCore& dindel)
        : GermlineIndelAlleleInfo(indelKey, indelData, indelReportInfo)
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


std::ostream& operator<<(std::ostream& os,const GermlineDiploidSiteAlleleInfo& smod);

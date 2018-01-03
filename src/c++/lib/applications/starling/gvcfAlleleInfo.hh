//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2018 Illumina, Inc.
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
    unsigned lead,
    unsigned trail,
    const unsigned commonPrefixLength,
    const IndelKey& indelKey,
    ALIGNPATH::path_t& cigar);




/// this object contains information shared by any germline variant allele
///
/// variant here means SNV or indel
///
/// SimpleGenotype here means information pertaining to the core genotyping algorithm for SNVs or Indels (where for indels this genotype is
///    restricted to a single alt allele). This information is supplemented with additional details required for full VCF records in the "Call"
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
        isForcedOutput = false;
        strandBias = 0.;
    }

    /// if true, this allele must appear in the VCF output:
    bool isForcedOutput = false;

    /// does this allele appear biased to one strand over all samples?:
    double strandBias = 0.;
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


/// restrict to the case where variant is site/SNV
struct GermlineSiteAlleleInfo : public GermlineVariantAlleleInfo
{
    typedef GermlineVariantAlleleInfo base_t;

    explicit
    GermlineSiteAlleleInfo(const BASE_ID::index_t initBaseId = BASE_ID::ANY)
        : baseIndex(initBaseId)
    {}

    void
    clear()
    {
        base_t::clear();
        baseIndex = BASE_ID::ANY;
    }

    BASE_ID::index_t baseIndex = BASE_ID::ANY;
};

std::ostream& operator<<(std::ostream& os,const GermlineSiteAlleleInfo& allele);

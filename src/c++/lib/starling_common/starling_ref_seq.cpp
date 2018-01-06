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

/// \file
/// \author Chris Saunders
///

#include "starling_common/starling_ref_seq.hh"

#include "common/Exceptions.hh"
#include "htsapi/samtools_fasta_util.hh"
#include "htsapi/bam_header_util.hh"



void
setRefSegment(
    const starling_base_options& opt,
    const std::string& chrom,
    const known_pos_range2& range,
    reference_contig_segment& ref)
{
    assert(! chrom.empty());

    ref.set_offset(range.begin_pos());
    // note: the ref function below takes closed-closed endpoints, so we subtract one from endPos
    get_standardized_region_seq(opt.referenceFilename, chrom, range.begin_pos(), range.end_pos()-1, ref.seq());
}



std::string
getSamtoolsRegionString(
    const std::string& chromName,
    const known_pos_range2& range)
{
    std::ostringstream oss;
    oss << chromName << ':' << range.begin_pos()+1 << '-' << range.end_pos();
    return oss.str();
}



void
getStrelkaAnalysisRegionInfo(
    const std::string& chrom,
    const int32_t beginPos,
    const int32_t endPos,
    const unsigned supplementalRegionBorderSize,
    AnalysisRegionInfo& rinfo)
{
    rinfo.regionChrom = chrom;
    rinfo.regionRange.set_range(beginPos, endPos);

    static const pos_t fixedRegionBorderSize(100);
    rinfo.streamerRegionRange = rinfo.regionRange;
    rinfo.streamerRegionRange.expandBy(fixedRegionBorderSize + supplementalRegionBorderSize);
    rinfo.streamerRegionRange.makeNonNegative();

    rinfo.streamerRegion = getSamtoolsRegionString(rinfo.regionChrom, rinfo.streamerRegionRange);

    static const pos_t refRegionPadSize(400);
    rinfo.refRegionRange = rinfo.streamerRegionRange;
    rinfo.refRegionRange.expandBy(refRegionPadSize);
    rinfo.refRegionRange.makeNonNegative();
}



void
getStrelkaAnalysisRegionInfo(
    const std::string& region,
    const unsigned supplementalRegionBorderSize,
    AnalysisRegionInfo& rinfo)
{
    int32_t regionBeginPos(0), regionEndPos(0);
    parse_bam_region(region.c_str(), rinfo.regionChrom, regionBeginPos, regionEndPos);
    getStrelkaAnalysisRegionInfo(rinfo.regionChrom, regionBeginPos, regionEndPos, supplementalRegionBorderSize, rinfo);
}



void
getStrelkaAnalysisRegions(
    const starling_base_options& opt,
    const std::string& referenceAlignmentFilename,
    const bam_header_info& referenceHeaderInfo,
    const unsigned supplementalRegionBorderSize,
    std::vector<AnalysisRegionInfo>& regionInfoList)
{
    /// TODO add test that regions do not intersect
    const unsigned regionCount(opt.regions.size());
    regionInfoList.resize(regionCount);

    for (unsigned regionIndex(0); regionIndex < regionCount; ++regionIndex)
    {
        const std::string& region(opt.regions[regionIndex]);
        auto& regionInfo(regionInfoList[regionIndex]);
        getStrelkaAnalysisRegionInfo(region, supplementalRegionBorderSize, regionInfo);

        // check that target region chrom exists in bam headers:
        if (not referenceHeaderInfo.chrom_to_index.count(regionInfo.regionChrom))
        {
            using namespace illumina::common;
            std::ostringstream oss;
            oss << "Region contig name: '" << regionInfo.regionChrom
                << "' is not found in the header of BAM/CRAM file: '" << referenceAlignmentFilename
                << "'";
            BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
        }
    }
}

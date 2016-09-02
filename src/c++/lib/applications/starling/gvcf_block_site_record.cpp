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

#include "gvcf_block_site_record.hh"
#include "blt_util/compat_util.hh"



static
bool
check_block_single_tolerance(const stream_stat& ss,
                             const int min,
                             const int tol)
{
    return ((min + tol) >= ss.max()/2.0);       // hack to get nova vcfs into acceptable size ramge, make check less stringent across the board
}



static
bool
check_block_tolerance(const stream_stat& ss,
                      const double frac_tol,
                      const int abs_tol)
{
    const int min(static_cast<int>(compat_round(ss.min())));
//    log_os << min << "\n";
//    return true;
    if (check_block_single_tolerance(ss,min,abs_tol)) return true;
    const int ftol(static_cast<int>(std::floor(min * frac_tol)));
    if (ftol <= abs_tol) return false;
    return check_block_single_tolerance(ss, min, ftol);
}



static
bool
is_new_value_blockable(const int new_val,
                       const stream_stat& ss,
                       const double frac_tol,
                       const int abs_tol,
                       const bool is_new_val = true,
                       const bool is_old_val = true)
{
    if (!(is_new_val && is_old_val)) return (is_new_val == is_old_val);

    stream_stat ss2(ss);
    ss2.add(new_val);
    return check_block_tolerance(ss2,frac_tol,abs_tol);
}

// sites that are 0/1 can be compressed if their non-ref allele ratios are low enough. So, fix those here
static const char* map_gt_to_homref(const char* gt)
{
    if (0 == std::strcmp("0/1", gt))
        return "0/0";
    return gt;
}



bool
gvcf_block_site_record::
testCanSiteJoinSampleBlockShared(
    const GermlineSiteLocusInfo& locus,
    const unsigned sampleIndex) const
{
    // pos must be +1 from end of record:
    if ((pos+count) != locus.pos) return false;

    const LocusSampleInfo& sampleInfo(getSample(0));
    const LocusSampleInfo& inputSampleInfo(locus.getSample(sampleIndex));

    // filters must match:
    if (not (filters == locus.filters)) return false;
    if (not (sampleInfo.filters == inputSampleInfo.filters)) return false;

    if (is_nonref() || locus.is_nonref()) return false;

    if (! is_new_value_blockable(locus.n_used_calls,
                                 block_dpu,frac_tol,abs_tol))
    {
        return false;
    }
    if (! is_new_value_blockable(locus.n_unused_calls,
                                 block_dpf,frac_tol,abs_tol))
    {
        return false;
    }

    return true;
}



void
gvcf_block_site_record::
joinSiteToSampleBlockShared(
    const GermlineSiteLocusInfo& locus,
    const unsigned /*sampleIndex*/)
{
    //LocusSampleInfo& sampleInfo(getSample(sampleIndex));
    //const LocusSampleInfo& inputSampleInfo(si.getSample(sampleIndex));

    if (count == 0)
    {
        pos = locus.pos;
        ref = locus.ref;
    }

    block_dpu.add(locus.n_used_calls);
    block_dpf.add(locus.n_unused_calls);

    count += 1;
}



bool
gvcf_block_site_record::
testCanSiteJoinSampleBlock(
    const GermlineDiploidSiteLocusInfo& locus,
    const unsigned sampleIndex) const
{
    if (count==0) return true;

    if (not testCanSiteJoinSampleBlockShared(locus,sampleIndex)) return false;

    const LocusSampleInfo& inputSampleInfo(locus.getSample(sampleIndex));

    if (gt != map_gt_to_homref(locus.get_gt())) return false;

    // coverage states must match:
    if (is_covered != locus.allele.is_covered) return false;
    if (is_used_covered != locus.allele.is_used_covered) return false;

    // ploidy must match
    if (ploidy != locus.dgt.ploidy) return false;

    // test blocking values:
    bool is_gqx(not locus.isRefUnknown());
    if (is_gqx)
    {
        is_gqx = locus.allele.is_gqx_tmp();
    }

    if (! is_new_value_blockable(inputSampleInfo.gqx,
                                 block_gqx,frac_tol,abs_tol,
                                 is_gqx,
                                 has_call))
    {
        return false;
    }

    return true;
}



void
gvcf_block_site_record::
joinSiteToSampleBlock(
    const GermlineDiploidSiteLocusInfo& locus,
    const unsigned sampleIndex)
{
    LocusSampleInfo& sampleInfo(getSample(0));
    const LocusSampleInfo& inputSampleInfo(locus.getSample(sampleIndex));

    bool is_gqx(not locus.isRefUnknown());
    if (is_gqx)
    {
        is_gqx = locus.allele.is_gqx_tmp();
    }

    if (count == 0)
    {
        filters = locus.filters;
        sampleInfo.filters = inputSampleInfo.filters;
        setNonRef(locus.is_nonref());
        gt = map_gt_to_homref(locus.get_gt());
        is_used_covered = locus.allele.is_used_covered;
        is_covered = locus.allele.is_covered;
        ploidy = locus.dgt.ploidy;
        has_call = is_gqx;
    }

    if (is_gqx)
    {
        block_gqx.add(inputSampleInfo.gqx);
    }

    joinSiteToSampleBlockShared(locus, sampleIndex);
}



bool
gvcf_block_site_record::
testCanSiteJoinSampleBlock(
    const GermlineContinuousSiteLocusInfo& locus,
    const unsigned sampleIndex) const
{
    if (locus.altAlleles.size() != 1)
        return false;

    if (count==0) return true;

    if (has_call && locus.altAlleles.empty())
        return false;

    // ploidy must match. This catches mixing digt && continuous
    if (ploidy != -1) return false;

    if(not testCanSiteJoinSampleBlockShared(locus,sampleIndex)) return false;

    const LocusSampleInfo& inputSampleInfo(locus.getSample(sampleIndex));

    if (gt != map_gt_to_homref(locus.get_gt(locus.altAlleles.front()))) return false;

    // coverage states must match:
    if (is_covered != (locus.n_used_calls != 0 || locus.n_unused_calls != 0)) return false;
    if (is_used_covered != (locus.n_used_calls != 0)) return false;

    if (has_call)
    {
        // test blocking values:
        if (! is_new_value_blockable(inputSampleInfo.gqx,
                                     block_gqx,frac_tol,abs_tol,
                                     true,
                                     true))
        {
            return false;
        }
    }

    return true;
}



void
gvcf_block_site_record::
joinSiteToSampleBlock(
    const GermlineContinuousSiteLocusInfo& locus,
    const unsigned sampleIndex)
{
    LocusSampleInfo& sampleInfo(getSample(0));
    const LocusSampleInfo& inputSampleInfo(locus.getSample(sampleIndex));

    if (count == 0)
    {
        if (!locus.altAlleles.empty())
        {
            filters = locus.filters;
            sampleInfo.filters = inputSampleInfo.filters;
            gt = map_gt_to_homref(locus.get_gt(locus.altAlleles.front()));
            setNonRef(locus.is_nonref());
            // TODO: handle no coverage regions in continuous
            has_call = true;
        }
        else
        {
            has_call = false;
        }
        is_used_covered = locus.n_used_calls != 0;
        is_covered = locus.n_used_calls != 0 || locus.n_unused_calls != 0;
        ploidy = -1;
    }

    // TODO: handle no coverage regions in continuous
    if (locus.altAlleles.size() == 1)
    {
        block_gqx.add(inputSampleInfo.gqx);
    }

    joinSiteToSampleBlockShared(locus, sampleIndex);
}

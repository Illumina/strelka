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


#include "gvcf_locus_info.hh"

#include <iostream>
#include <map>
#include <typeinfo>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions.hpp>

void
shared_call_info::
write_filters(std::ostream& os) const
{
    if (filters.none())
    {
        os << "PASS";
        return;
    }

    bool is_sep(false);
    for (unsigned i(0); i<VCF_FILTERS::SIZE; ++i)
    {
        if (! filters.test(i)) continue;

        if (is_sep)
        {
            os << ";";
        }
        else
        {
            is_sep=true;
        }
        os << VCF_FILTERS::get_label(i);
    }
}


pos_t digt_indel_info::end() const
{
    pos_t result = 0;
    for (auto& x : _calls)
        result = std::max(result, x._ik.right_pos());
    return result;
}

void shared_indel_call_info::set_hap_cigar(
    const unsigned lead,
    const unsigned trail)
{
    using namespace ALIGNPATH;

    cigar.clear();
    if (lead)
    {
        cigar.push_back(path_segment(MATCH,lead));
    }
    if (_ik.delete_length())
    {
        cigar.push_back(path_segment(DELETE,_ik.delete_length()));
    }
    if (_ik.insert_length())
    {
        cigar.push_back(path_segment(INSERT,_ik.insert_length()));
    }
    if (trail)
    {
        cigar.push_back(path_segment(MATCH,trail));
    }
}


std::map<std::string, double>
digt_indel_info::get_indel_qscore_features(
    const double chrom_depth) const
{
    const auto& call(first());
    std::map<std::string, double> res;
    res["QUAL"]             = call._dindel.indel_qphred /(1.*chrom_depth);
    res["F_GQX"]            = call.gqx /(1.*chrom_depth);
    res["F_GQ"]             = call.gq /(1.*chrom_depth); // N.B. Not used at time of writing; normalization uncertain
    res["REFREP1"]          = call._iri.ref_repeat_count;

    res["IDREP1"]           = call._iri.indel_repeat_count;
    res["RULEN1"]           = call._iri.repeat_unit.length(); //isri.depth;               //This feature actually means the length of the RU string

    unsigned ref_count(0);
    ref_count = std::max(ref_count,first()._isri.n_q30_ref_reads);

    const double r0 = ref_count;
    const double r1 = call._isri.n_q30_indel_reads;
    const double r2 = call._isri.n_q30_alt_reads;
    res["AD0"]              = r0/(1.0*chrom_depth);
    res["AD1"]              = r1/(1.0*chrom_depth);
    res["AD2"]              = r2/(1.0*chrom_depth);
    // allele bias metrics
    // cdf of binomial prob of seeing no more than the number of 'allele A' reads out of A reads + B reads, given p=0.5
    double allelebiaslower = cdf(boost::math::binomial(r0+r1,0.5),r0);
    // cdf of binomial prob of seeing no more than the number of 'allele B' reads out of A reads + B reads, given p=0.5
    double allelebiasupper = cdf(boost::math::binomial(r0+r1,0.5),r1);
    if ( is_hetalt() )
    {
        allelebiaslower = cdf(boost::math::binomial(r2+r1,0.5),r1);
        allelebiasupper = cdf(boost::math::binomial(r2+r1,0.5),r2);
    }
    res["ABlower"]          = -std::log(allelebiaslower+1.e-30); // +1e-30 to avoid log(0) in extreme cases
    res["AB"]               = -std::log(std::min(1.,2.*std::min(allelebiaslower,allelebiasupper))+1.e-30);

    res["F_DPI"]            = call._isri.depth/(1.0*chrom_depth);
    return res;
}


static void add_cigar_to_ploidy(
    const ALIGNPATH::path_t& apath,
    std::vector<unsigned>& ploidy)
{
    using namespace ALIGNPATH;
    int offset(-1);
    for (const auto& ps : apath)
    {
        if (is_segment_align_match(ps.type))
        {
            for (unsigned j(0); j<ps.length; ++j)
            {
                if (offset>=0) ploidy[offset]++;
                offset++;
            }
        }
        else if (ps.type==DELETE)
        {
            offset+=ps.length;
        }
    }
}


void digt_indel_info::add_overlap(const reference_contig_segment& ref, digt_indel_info& overlap)
{
    assert(_calls.size() == 1);
    assert(overlap._calls.size() == 1);

    auto& call(first());
    auto& overlap_call(overlap.first());



    // there's going to be 1 (possibly empty) fill range in front of one haplotype
    // and one possibly empty fill range on the back of one haplotype
    std::string leading_seq,trailing_seq;
    auto indel_end_pos=std::max(overlap_call._ik.right_pos(),call._ik.right_pos());

    const pos_t indel_begin_pos(pos-1);

    // add shared information (to the first indel only)
    // TODO: Reevaluate this. Since we're merging, we can do this on output
    // make extended vcf ref seq:
    std::string tmp;
    ref.get_substring(indel_begin_pos,(indel_end_pos-indel_begin_pos),tmp);
    call._iri.vcf_ref_seq = tmp;

    ploidy.resize(indel_end_pos-pos,0);

    auto munge_indel = [&] (digt_indel_info& ii)
    {
        auto& this_call(ii.first());
        // extend leading sequence start back 1 for vcf compat, and end back 1 to concat with vcf_indel_seq
        ref.get_substring(indel_begin_pos,(ii.pos-indel_begin_pos)-1,leading_seq);
        const unsigned trail_len(indel_end_pos-this_call._ik.right_pos());
        ref.get_substring(indel_end_pos-trail_len,trail_len,trailing_seq);


        this_call._iri.vcf_indel_seq = leading_seq + this_call._iri.vcf_indel_seq + trailing_seq;

        this_call.set_hap_cigar(leading_seq.size()+1,
                                trailing_seq.size());

        // add to the ploidy object:
        add_cigar_to_ploidy(this_call.cigar,this->ploidy);
    };
    munge_indel(*this);
    munge_indel(overlap);


    //reduce qual and gt to the lowest of the set:
    call._dindel.indel_qphred = std::min(call._dindel.indel_qphred, overlap.first()._dindel.indel_qphred);
    call._dindel.max_gt_qphred = std::min(call._dindel.max_gt_qphred, overlap.first()._dindel.max_gt_qphred);


    // combine filter flags
    call.filters |= overlap.first().filters;
    // combine QScores. Since the "unset" value is -1, this complex logic is necessary
    if (call.Qscore <0)
        call.Qscore = overlap.first().Qscore;
    else if (overlap.first().Qscore >= 0)
        call.Qscore = std::min(call.Qscore, overlap.first().Qscore);
    call.gqx = std::min(call.gqx, overlap.first().gqx);
    call.gq = std::min(call.gq, overlap.first().gq);

    _calls.push_back(overlap.first());
}




std::map<std::string, double>
digt_site_info::
get_site_qscore_features(
    double chrom_depth) const
{
    std::map<std::string, double> res;

    res["QUAL"]               = dgt.genome.snp_qphred / (1.*chrom_depth);
    res["F_GQX"]              = smod.gqx / (1.*chrom_depth);
    res["F_GQ"]               = smod.gq / (1.*chrom_depth);
    res["I_SNVSB"]            = smod.strand_bias;
    res["I_SNVHPOL"]          = hpol;

    //we need to handle the scaling of DP better for high depth cases
    res["F_DP"]               = n_used_calls/(1.0*chrom_depth);
    res["F_DPF"]              = n_unused_calls/(1.0*chrom_depth);
    res["AD0"]                = known_counts[dgt.ref_gt]/(1.0*chrom_depth);
    res["AD1"]                = 0.0;          // set below

    res["I_MQ"]               = MQ;
    res["I_ReadPosRankSum"]   = ReadPosRankSum;
    res["I_BaseQRankSum"]     = BaseQRankSum;
    res["I_MQRankSum"]        = MQRankSum;
    res["I_RawPos"]           = rawPos;         //the average position value within a read of alt allele
    res["I_RawBaseQ"]         = avgBaseQ;       //The average baseQ of the position of alt allele
    for (unsigned b(0); b<N_BASE; ++b)
    {
        if (b==dgt.ref_gt) continue;
        if (DIGT::expect2(b,smod.max_gt))
        {
            res["AD1"] =  known_counts[b]/(1.0*chrom_depth);
            // allele bias metrics
            double r0 = known_counts[dgt.ref_gt];
            double r1 = known_counts[b];
            double allelebiaslower  = cdf(boost::math::binomial(r0+r1,0.5),r0);
            double allelebiasupper  = cdf(boost::math::binomial(r0+r1,0.5),r1);
            res["ABlower"]          = -log(allelebiaslower+1.e-30); // +1e-30 to avoid log(0) in extreme cases
            res["AB"]               = -log(std::min(1.,2.*std::min(allelebiaslower,allelebiasupper))+1.e-30);
        }
    }
    if ((res["F_DP"]+res["F_DPF"])>0.0)
    {
        res["VFStar"]           = res["AD1"]/(res["DP"]+res["DPF"]); //VFStar = AD2/(DP+DPF);
    }
    else
    {
        res["VFStar"]           = res["AD1"]/(1.0*chrom_depth); //default hack for
    }
    return res;
}






std::ostream&
operator<<(std::ostream& os,
           const shared_call_info& shmod)
{
    os << "gqx: " << shmod.gqx
       << " gq: " << shmod.gq;
    if (typeid(shmod) == typeid(digt_indel_call))
    {
        auto imod = dynamic_cast<const digt_indel_call&>(shmod);

        os << " max_gt: " << DIGT::label(imod.max_gt);
    }

    os << " filters: ";
    shmod.write_filters(os);

    return os;
}

std::ostream&
operator<<(std::ostream& os,
           const digt_call_info& smod)
{
    os << static_cast<shared_call_info>(smod) << '\n';

    os << "is_unknown: " << smod.is_unknown;
    os << " is_covered: " << smod.is_covered;
    os << " is_used_coverage: " << smod.is_used_covered;
    os << " is_zero_ploidy: " << smod.is_zero_ploidy;

    if (smod.modified_gt != MODIFIED_SITE_GT::NONE)
    {
        os << " modgt: " << MODIFIED_SITE_GT::get_label(smod.modified_gt);
    }

    return os;
}

std::ostream&
operator<<(std::ostream& os,
           const digt_site_info& si)
{
    os << "pos: " << (si.pos+1) << " " << si.get_gt();
    return os;
}

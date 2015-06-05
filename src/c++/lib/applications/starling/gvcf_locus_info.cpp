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


#include "gvcf_locus_info.hh"

#include <iostream>
#include <map>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions.hpp>

void
shared_modifiers::
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


std::map<std::string, double>
indel_info::
get_indel_qscore_features(
    const double chrom_depth) const
{
    std::map<std::string, double> res;
    res["QUAL"]             = dindel.indel_qphred /(1.*chrom_depth);
    res["F_GQX"]            = imod().gqx /(1.*chrom_depth);
    res["F_GQ"]             = imod().gq /(1.*chrom_depth); // N.B. Not used at time of writing; normalization uncertain
    res["REFREP1"]          = iri().ref_repeat_count;

    res["IDREP1"]           = iri().indel_repeat_count;
    res["RULEN1"]           = iri().repeat_unit.length(); //isri.depth;               //This feature actually means the length of the RU string

    unsigned ref_count(0);
    ref_count = std::max(ref_count,isri().n_q30_ref_reads);

    const double r0 = ref_count;
    const double r1 = isri().n_q30_indel_reads;
    const double r2 = isri().n_q30_alt_reads;
    res["AD0"]              = r0/(1.0*chrom_depth);
    res["AD1"]              = r1/(1.0*chrom_depth);
    res["AD2"]              = r2/(1.0*chrom_depth);
    // allele bias metrics
    // cdf of binomial prob of seeing no more than the number of 'allele A' reads out of A reads + B reads, given p=0.5
    double allelebiaslower = cdf(boost::math::binomial(r0+r1,0.5),r0);
    // cdf of binomial prob of seeing no more than the number of 'allele B' reads out of A reads + B reads, given p=0.5
    double allelebiasupper = cdf(boost::math::binomial(r0+r1,0.5),r1);
    if ( _imod.size() > 1 )
    {
        allelebiaslower = cdf(boost::math::binomial(r2+r1,0.5),r1);
        allelebiasupper = cdf(boost::math::binomial(r2+r1,0.5),r2);
    }
    res["ABlower"]          = -std::log(allelebiaslower+1.e-30); // +1e-30 to avoid log(0) in extreme cases
    res["AB"]               = -std::log(std::min(1.,2.*std::min(allelebiaslower,allelebiasupper))+1.e-30);

    res["F_DPI"]            = isri().depth/(1.0*chrom_depth);
    return res;
}

void indel_info::set_hap_cigar(
        const unsigned lead,
        const unsigned trail)
{
    using namespace ALIGNPATH;
    ALIGNPATH::path_t& apath(_imod.front().cigar);

    apath.clear();
    if (lead)
    {
        apath.push_back(path_segment(MATCH,lead));
    }
    if (ik().delete_length())
    {
        apath.push_back(path_segment(DELETE,ik().delete_length()));
    }
    if (ik().insert_length())
    {
        apath.push_back(path_segment(INSERT,ik().insert_length()));
    }
    if (trail)
    {
        apath.push_back(path_segment(MATCH,trail));
    }
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

pos_t indel_info::end() const
{
    pos_t result = 0;
    for (auto x : _ik)
        result = std::max(result, x.right_pos());
    return result;
}

void indel_info::add_overlap(const reference_contig_segment& ref, indel_info& overlap)
{
    assert(_imod.size() == 1);

    // there's going to be 1 (possibly empty) fill range in front of one haplotype
    // and one possibly empty fill range on the back of one haplotype
    std::string leading_seq,trailing_seq;
    auto indel_end_pos=std::max(overlap.ik().right_pos(),ik().right_pos());

    const pos_t indel_begin_pos(pos-1);

    // add shared information (to the first indel only)
    // TODO: Reevaluate this. Since we're merging, we can do this on output
    // make extended vcf ref seq:
    std::string tmp;
    ref.get_substring(indel_begin_pos,(indel_end_pos-indel_begin_pos),tmp);
    _iri.front().vcf_ref_seq = tmp;

    imod().ploidy.resize(indel_end_pos-pos,0);

    auto munge_indel = [&] (indel_info& ii)
    {
        // extend leading sequence start back 1 for vcf compat, and end back 1 to concat with vcf_indel_seq
        ref.get_substring(indel_begin_pos,(ii.pos-indel_begin_pos)-1,leading_seq);
        const unsigned trail_len(indel_end_pos-ii.ik().right_pos());
        ref.get_substring(indel_end_pos-trail_len,trail_len,trailing_seq);


        ii._iri.front().vcf_indel_seq = leading_seq + ii.iri().vcf_indel_seq + trailing_seq;

        ii.set_hap_cigar(leading_seq.size()+1,
                      trailing_seq.size());

        // add to the ploidy object:
        add_cigar_to_ploidy(ii.imod().cigar,this->_imod.front().ploidy);
    };
    munge_indel(*this);
    munge_indel(overlap);


    //reduce qual and gt to the lowest of the set:
    this->dindel.indel_qphred = std::min(this->dindel.indel_qphred, overlap.dindel.indel_qphred);
    this->dindel.max_gt_qphred = std::min(this->dindel.max_gt_qphred, overlap.dindel.max_gt_qphred);


    // combine filter flags
    this->imod().filters |= overlap.imod().filters;
    // combine QScores. Since the "unset" value is -1, this complex logic is necessary
    if (this->imod().Qscore <0)
        this->imod().Qscore = overlap.imod().Qscore;
    else if (overlap.imod().Qscore >= 0)
        this->imod().Qscore = std::min(this->imod().Qscore, overlap.imod().Qscore);

    _ik.push_back(overlap.ik());
    _iri.push_back(overlap.iri());
    _isri.push_back(overlap.isri());
    _imod.push_back(overlap.imod());
}




std::map<std::string, double>
site_info::
get_site_qscore_features(
    double chrom_depth) const
{
    std::map<std::string, double> res;

    res["QUAL"]               = dgt.genome.snp_qphred / (1.*chrom_depth);
    res["F_GQX"]              = smod.gqx / (1.*chrom_depth);
    res["F_GQ"]               = smod.gq / (1.*chrom_depth);
    res["I_SNVSB"]            = dgt.sb;
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
           const shared_modifiers& shmod)
{
    os << "gqx: " << shmod.gqx
       << " gq: " << shmod.gq
       << " max_gt: " << DIGT::label(shmod.max_gt);

    os << " filters: ";
    shmod.write_filters(os);

    return os;
}

std::ostream&
operator<<(std::ostream& os,
           const site_modifiers& smod)
{
    os << static_cast<shared_modifiers>(smod) << '\n';

    os << "is_unknown: " << smod.is_unknown;
    os << " is_covered: " << smod.is_covered;
    os << " is_used_coverage: " << smod.is_used_covered;
    os << " is_zero_ploidy: " << smod.is_zero_ploidy;
    os << " is_block: " << smod.is_block;

    if (smod.modified_gt != MODIFIED_SITE_GT::NONE)
    {
        os << " modgt: " << MODIFIED_SITE_GT::get_label(smod.modified_gt);
    }

    return os;
}

std::ostream&
operator<<(std::ostream& os,
           const site_info& si)
{
    os << "pos: " << (si.pos+1) << " " << si.get_gt();
    return os;
}

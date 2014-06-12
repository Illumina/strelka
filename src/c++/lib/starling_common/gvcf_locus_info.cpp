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

/// \file
///
/// \author Chris Saunders
///


#include "starling_common/gvcf_locus_info.hh"

#include <iostream>
#include <map>


void
shared_modifiers::
write_filters(std::ostream& os) const {

    if (filters.none()) {
        os << "PASS";
        return;
    }

    bool is_sep(false);
    for (unsigned i(0); i<VCF_FILTERS::SIZE; ++i) {
        if (filters.test(i)) {
            if (is_sep) { os << ";"; }
            else       { is_sep=true; }
            os << VCF_FILTERS::get_label(i);
        }
    }
}


std::map<std::string, double> indel_info::get_qscore_features() {
    this->calc_vqsr_metrics();

    // set GQ and GQX
    if (dindel.max_gt != dindel.max_gt_poly) {
        imod.gqx=0;
    } else {
        imod.gqx=std::min(dindel.max_gt_poly_qphred,dindel.max_gt_qphred);
    }
    imod.max_gt=dindel.max_gt_poly;
    imod.gq=dindel.max_gt_poly_qphred;
    //F_DPI   F_GQ    F_GQX   CLASS   RULEN1  REFREP1 IDREP1  AD0     AD1

    std::map<std::string, double> res;
    res["QUAL"]             = dindel.indel_qphred;
    res["F_GQX"]            = imod.gqx;
    res["F_GQ"]             = imod.gq;
    res["REFREP1"]          = iri.ref_repeat_count;

    //res["LENGTH"]           = ik.length;
    res["IDREP1"]           = iri.indel_repeat_count;
    res["RULEN1"]           = iri.repeat_unit.length(); //isri.depth;               //This feature actually means the length of the RU string

    if (imod.is_overlap) {
        // hack for overlap case
        //res["REFREP2"]          = iri.ref_repeat_count;
        //res["IDREP2"]           = iri.indel_repeat_count;
        //res["RULEN2"]           = iri.repeat_unit.length();
    }
    unsigned ref_count(0);
    ref_count = std::max(ref_count,isri.n_q30_ref_reads);
    res["AD0"]              = ref_count;
    res["AD1"]              = isri.n_q30_indel_reads;
    res["F_DPI"]            = isri.depth;
//    res["MQ"]               = MQ;
//    res["ReadPosRankSum"]   = ReadPosRankSum;
//    res["BaseQRankSum"]     = BaseQRankSum;
//    res["MQRankSum"]        = MQRankSum;
    return res;
}

void indel_info::calc_vqsr_metrics() {
    this->MQ                = 0.0; //this->ik.mapq_val*1.0/this->ik.mapq_n;
    this->ReadPosRankSum    = 1.0;
    this->MQRankSum         = 2.0;
    this->BaseQRankSum      = 3.0;
//    this->re              = this->ik.mapq_val*1.0/this->ik.mapq_n;
}


std::map<std::string, double> site_info::get_qscore_features() {
    std::map<std::string, double> res;
    res["QUAL"]               = dgt.genome.snp_qphred;;
    res["F_GQX"]              = smod.gqx;
    res["F_GQ"]               = smod.gq;
    res["I_SNVSB"]            = dgt.sb;
    res["I_SNVHPOL"]          = hpol;
    res["F_DP"]               = n_used_calls;
    res["F_DPF"]              = n_unused_calls;
    res["AD0"]                = known_counts[dgt.ref_gt];
    res["AD1"]                = 0.0;          // set below
    res["I_MQ"]               = MQ;
    res["I_ReadPosRankSum"]   = ReadPosRankSum;
    res["I_BaseQRankSum"]     = BaseQRankSum;
    res["I_MQRankSum"]        = MQRankSum;
    for (unsigned b(0); b<N_BASE; ++b) {
        if (b==dgt.ref_gt) continue;
        if (DIGT::expect2(b,smod.max_gt))
            res["AD1"] =  known_counts[b];
    }
    if ((res["F_DP"]+res["F_DPF"])>0.0) {
        res["VFStar"]           = res["AD1"]/(res["DP"]+res["DPF"]); //VFStar = AD2/(DP+DPF);
    }
    else {
        res["VFStar"]           = res["AD1"]/(30.0); //default hack for
    }
    return res;
}


std::ostream&
operator<<(std::ostream& os,
           const shared_modifiers& shmod) {

    os << "gqx: " << shmod.gqx
       << " gq: " << shmod.gq
       << " max_gt: " << DIGT::label(shmod.max_gt);

    os << " filters: ";
    shmod.write_filters(os);

    return os;
}

std::ostream&
operator<<(std::ostream& os,
           const site_modifiers& smod) {

    os << static_cast<shared_modifiers>(smod) << '\n';

    os << "is_unknown: " << smod.is_unknown;
    os << " is_covered: " << smod.is_covered;
    os << " is_used_coverage: " << smod.is_used_covered;
    os << " is_zero_ploidy: " << smod.is_zero_ploidy;
    os << " is_block: " << smod.is_block;

    if (smod.modified_gt != MODIFIED_SITE_GT::NONE) {
        os << " modgt: " << MODIFIED_SITE_GT::get_label(smod.modified_gt);
    }

    return os;
}

std::ostream&
operator<<(std::ostream& os,
           const site_info& si) {
    os << "pos: " << (si.pos+1) << " " << si.get_gt();
    return os;
}

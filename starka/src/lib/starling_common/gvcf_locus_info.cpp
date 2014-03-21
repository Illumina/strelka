// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2013 Illumina, Inc.
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
    std::map<std::string, double> res;
    res["GQX"]              = imod.gqx;
    res["GQ"]               = imod.gq;
    res["REFREP"]           = iri.ref_repeat_count;
    res["LENGTH"]           = ik.length;
    res["IDREP"]            = iri.indel_repeat_count;
    res["RU"]               = iri.repeat_unit.length(); //isri.depth;               //This feature actually means the length of the RU string
    res["MQ"]               = MQ;
    res["ReadPosRankSum"]   = ReadPosRankSum;
    res["BaseQRankSum"]     = BaseQRankSum;
    res["MQRankSum"]        = MQRankSum;
    return res;
}

void indel_info::calc_vqsr_metrics(){
    this->MQ                = 0.0; //this->ik.mapq_val*1.0/this->ik.mapq_n;
    this->ReadPosRankSum    = 1.0;
    this->MQRankSum         = 2.0;
    this->BaseQRankSum      = 3.0;
//    this->re              = this->ik.mapq_val*1.0/this->ik.mapq_n;
}


std::map<std::string, double> site_info::get_qscore_features() {
    std::map<std::string, double> res;
    res["GQX"]              = smod.gqx;
    res["GQ"]               = smod.gq;
    res["SNVSB"]            = dgt.sb;
    res["SNVHPOL"]          = hpol;
    res["DP"]               = n_used_calls;
    res["DPF"]              = n_unused_calls;
    res["AD"]               = known_counts[dgt.ref_gt];
    res["AD2"]              = 0.0;          // set below
    res["MQ"]               = MQ;
    res["ReadPosRankSum"]   = ReadPosRankSum;
    res["BaseQRankSum"]     = BaseQRankSum;
    res["MQRankSum"]        = MQRankSum;
    for (unsigned b(0); b<N_BASE; ++b) {
        if (b==dgt.ref_gt) continue;
        if (DIGT::expect2(b,smod.max_gt))
            res["AD2"] =  known_counts[b];
    }
    if ((res["DP"]+res["DPF"])>0.0) {
        res["VFStar"]           = res["AD2"]/(res["DP"]+res["DPF"]); //VFStar = AD2/(DP+DPF);
    }
    else {
        res["VFStar"]           = res["AD2"]/(30.0); //default hack for
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

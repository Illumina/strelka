// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2012 Illumina, Inc.
//
// This software is covered by the "Illumina Genome Analyzer Software
// License Agreement" and the "Illumina Source Code License Agreement",
// and certain third party copyright/licenses, and any user of this
// source file is bound by the terms therein (see accompanying files
// Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
// Illumina_Source_Code_License_Agreement.pdf and third party
// copyright/license notices).
//
//

/// \file
///
/// \author Chris Saunders
///

#include "gvcf_aggregator.hh"

#include <iostream>


#if 0
VcfRecord {
    reset(site_info);
    reset(site_info,site_modifiers);
    reset(indel_info);
};


gvcf_blocker {

    add(const site_info&,
        const site_modifiers&);

    add(const indel_info&, indel_modifiers&);

};
#endif



gvcf_aggregator::
gvcf_aggregator(const starling_options& opt,
                const pos_range& report_range,
                const reference_contig_segment& ref,
                std::ostream* osptr)
    : _opt(opt)
    , _report_range(report_range.begin_pos,report_range.end_pos)
    , _ref(ref)
    , _osptr(osptr)
    , _chrom(NULL)
    , _indel_end_pos(0)
    , _indel_buffer_size(0)
    , _site_buffer_size(0)
{
    assert(report_range.is_begin_pos);
    assert(report_range.is_end_pos);

    if(opt.is_gvcf_output()) {
        assert(NULL != osptr);
    }
}



gvcf_aggregator::
~gvcf_aggregator() {
    process_overlaps();
}



void
gvcf_aggregator::
add_site(const pos_t pos,
         const char ref,
         const unsigned n_used_calls,
         const unsigned n_unused_calls,
         const snp_pos_info& good_pi,
         const diploid_genotype& dgt,
         const bool is_nf_snp,
         const double sb,
         const unsigned hpol) {

    if(0 != _indel_buffer_size) {
        if(pos>=_indel_end_pos) {
            process_overlaps();
        } else {
            while(_site_buffer.size() <= _site_buffer_size) {
                _site_buffer.push_back(site_info());
            }
            _site_buffer[_site_buffer_size++].init(pos,ref,n_used_calls,n_unused_calls,good_pi,dgt,is_nf_snp,sb,hpol);
            return;
        }
    }

    // write_site:
    *_osptr << _chrom << "\t" << (pos+1) << "\t" << "SNP" << "\t" << ref << "\n";
}



static
bool
is_het_indel(const starling_diploid_indel_core& dindel) {
    return (dindel.max_gt==STAR_DIINDEL::HET);
}

static
bool
is_no_indel(const starling_diploid_indel_core& dindel) {
    return (dindel.max_gt==STAR_DIINDEL::NOINDEL);
}



void
gvcf_aggregator::
add_indel(const pos_t pos,
          const indel_key ik,
          const starling_diploid_indel_core& dindel,
          const starling_indel_report_info& iri,
          const starling_indel_sample_report_info& isri) {

    // we can't handle breakends at all right now:
    if(ik.is_breakpoint()) return;

    // don't handle max_gt=="ref" cases for now:
    if(is_no_indel(dindel)) return;

    // check to see if we add this indel to the buffer:
    if(0 != _indel_buffer_size) {
        // check if this indel overlaps the buffer -- note this deleberatly picks up adjacent deletions:
        if(pos<=_indel_end_pos) {
            _indel_end_pos=std::max(_indel_end_pos,ik.right_pos());
        } else {
            process_overlaps();
        }
    }

    while(_indel_buffer.size() <= _indel_buffer_size) {
        _indel_buffer.push_back(indel_info());
    }
    _indel_buffer[_indel_buffer_size++].init(pos,ik,dindel,iri,isri);
    _indel_end_pos=ik.right_pos();
}



static
bool
is_simple_indel_overlap(const std::vector<indel_info>& indel_buffer,
                        const unsigned size) {

    return (size==2 &&
            is_het_indel(indel_buffer[0].dindel) &&
            is_het_indel(indel_buffer[1].dindel));

}


static
void
get_hap_cigar(const std::string lead,
              const std::string ref,
              const std::string alt,
              const std::string trail,
              ALIGNPATH::path_t& apath) {

    using namespace ALIGNPATH;

    apath.clear();
    if(lead.size()) {
        apath.push_back(path_segment(MATCH,lead.size()));
    }
    if(ref.size()) {
        apath.push_back(path_segment(DELETE,ref.size()));
    }
    if(alt.size()) {
        apath.push_back(path_segment(INSERT,alt.size()));
    }
    if(trail.size()) {
        apath.push_back(path_segment(MATCH,trail.size()));
    }
}



static
void
get_hap_cigar(const std::string ref,
              const std::string alt,
              ALIGNPATH::path_t& apath) {

    static const std::string lead("N");
    static const std::string trail("");
    get_hap_cigar(lead,ref,alt,trail,apath);
}



// set the CIGAR string:
void
gvcf_aggregator::
modify_single_indel_record() {
    assert(_indel_buffer_size==1);

    indel_info& ii(_indel_buffer[0]);
    get_hap_cigar(ii.iri.ref_seq,ii.iri.indel_seq,ii.imod.cigar);
}



void
gvcf_aggregator::
modify_overlap_indel_record() {
    
    // can only handle simple 2-indel overlaps right now:
    assert(_indel_buffer_size==2);

    // accumutate all modification info in the *first* indel record:
    indel_info& ii(_indel_buffer[0]);

    ii.imod.is_overlap=true;

    // there's going to be 1 (possibly empty) fill range in front of one haplotype
    // and one possibly empty fill range on the back of one haplotype
    std::string leading_seq,trailing_seq;

    const pos_t indel_begin_pos(ii.pos-1);

    // add shared information (to the first indel only)
    // make extended vcf ref seq:
    _ref.get_substring(indel_begin_pos,(_indel_end_pos-indel_begin_pos),ii.iri.vcf_ref_seq);


    // add per-haplotype information:
    for(unsigned hap(0);hap<2;++hap) {
        //reduce qual and gt to the lowest of the set:
        if(hap) {
            if(ii.dindel.indel_qphred>_indel_buffer[hap].dindel.indel_qphred) {
                ii.dindel.indel_qphred = _indel_buffer[hap].dindel.indel_qphred;
            }
            if(ii.dindel.max_gt_qphred>_indel_buffer[hap].dindel.max_gt_qphred) {
                ii.dindel.max_gt_qphred = _indel_buffer[hap].dindel.max_gt_qphred;
            }
        }
        
        _ref.get_substring(indel_begin_pos,_indel_buffer[hap].pos-indel_begin_pos,leading_seq);
        const unsigned trail_len(_indel_end_pos-_indel_buffer[hap].ik.right_pos());
        _ref.get_substring(_indel_end_pos-trail_len,trail_len,trailing_seq);

        _indel_buffer[hap].iri.vcf_indel_seq = leading_seq + _indel_buffer[hap].iri.indel_seq + trailing_seq;

        get_hap_cigar(leading_seq,
                      _indel_buffer[hap].iri.ref_seq,
                      _indel_buffer[hap].iri.indel_seq,
                      trailing_seq,
                      _indel_buffer[hap].imod.cigar);
    }
}



static
void
add_indel_modifiers(const starling_options& opt,
                    indel_info& ii) {

    ii.imod.gqx=std::min(ii.dindel.indel_qphred,ii.dindel.max_gt_qphred);
    if(opt.is_gvcf_min_gqx) {
        if(ii.imod.gqx<opt.gvcf_min_gqx) ii.imod.set_filter(VCF_FILTERS::LowGQX);
    }

    if(opt.is_gvcf_max_depth) {
        if(ii.isri.depth > opt.gvcf_max_depth) ii.imod.set_filter(VCF_FILTERS::HighDepth);
    }
}



void
gvcf_aggregator::
write_indel_record() {

    assert(_indel_buffer_size>0);

    std::ostream& os(*_osptr);
    indel_info& ii(_indel_buffer[0]);

    // compute final mods:
    add_indel_modifiers(_opt,ii);

    os << _chrom << '\t'   // CHROM
       << ii.pos << '\t'   // POS
       << ".\t"            // ID
       << ii.iri.vcf_ref_seq << '\t'; // REF

    // ALT
    for(unsigned i(0);i<_indel_buffer_size;++i) {
        if(i) os << ',';
        os << ii.iri.vcf_indel_seq;
    }
    os << '\t';

    os << ii.dindel.indel_qphred << '\t'; //QUAL

    // FILTER:
    ii.imod.write_filters(os);
    os << '\t';

    // INFO
    os << "CIGAR=";
    for(unsigned i(0);i<_indel_buffer_size;++i) {
        if(i) os << ',';
        os << _indel_buffer[i].imod.cigar;
    }
    os << '\t';

    //FORMAT
    os << "GT:GQX" << '\t';

    //SAMPLE
    os << ii.get_gt() << ':' << ii.imod.gqx  << '\n';
}



void
gvcf_aggregator::
process_overlaps() {

    if(0==_indel_buffer_size) return;

    // do the overlap processing:
    if(_indel_buffer_size==1) {
        // simple case of no overlap:
        modify_single_indel_record();
    } else {
        if(is_simple_indel_overlap(_indel_buffer,_indel_buffer_size)){
            // handle the simplest possible overlap case (two hets):
            modify_overlap_indel_record();
        } else {
            // mark the whole region as conflicting
        }
    }

    write_indel_record();
    *_osptr << "INDEL_SIZE: " << _indel_buffer_size << "\n";

    // tmp debug:
    for(unsigned i(0);i<_site_buffer_size;++i) {
        *_osptr << _chrom << "\t" << (_site_buffer[i].pos+1) << "\t" << "BufferedSNP" << "\t" << _site_buffer[i].ref << "\n";
    }

    _indel_buffer_size = 0;
    _site_buffer_size = 0;
}




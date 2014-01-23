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

#include "blt_util/blt_exception.hh"
#include "starling_common/chrom_depth_map.hh"
#include "starling_common/gvcf_aggregator.hh"
#include "starling_common/gvcf_header.hh"


#include "boost/foreach.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>


#define DEBUG_GVCF


#ifdef DEBUG_GVCF
#include "blt_util/log.hh"
#endif



static
void
set_site_gt(const diploid_genotype::result_set& rs,
            site_modifiers& smod) {

    smod.max_gt=rs.max_gt;
    smod.gqx=rs.max_gt_qphred;
}

//legacy methof
static
void
set_site_filters(const gvcf_options& opt,
                 const gvcf_deriv_options& dopt,
                 site_info& si) {

      if (opt.is_min_gqx) {
            if (si.smod.gqx<opt.min_gqx) si.smod.set_filter(VCF_FILTERS::LowGQX);
        }

        if (dopt.is_max_depth) {
            if ((si.n_used_calls+si.n_unused_calls) > dopt.max_depth) si.smod.set_filter(VCF_FILTERS::HighDepth);
        }

        if (opt.is_max_base_filt) {
            const unsigned total_calls(si.n_used_calls+si.n_unused_calls);
            if (total_calls>0) {
                const double filt(static_cast<double>(si.n_unused_calls)/static_cast<double>(total_calls));
                if (filt>opt.max_base_filt) si.smod.set_filter(VCF_FILTERS::HighBaseFilt);
            }
        }

        if (si.dgt.is_snp) {
            if (opt.is_max_snv_sb) {
                if (si.dgt.sb>opt.max_snv_sb) si.smod.set_filter(VCF_FILTERS::HighSNVSB);
            }

            if (opt.is_max_snv_hpol) {
                if (static_cast<int>(si.hpol)>opt.max_snv_hpol) si.smod.set_filter(VCF_FILTERS::HighSNVHPOL);
            }
        }
}

void
set_site_filters_CM(const gvcf_options& opt,
                 const gvcf_deriv_options& dopt,
                 site_info& si,
                 calibration_models& model) {
    // Code for old command-line parameterized filter behaviour has been moved to calibration_models.cpp
    model.clasify_site(opt,dopt,si);
}



static
void
add_site_modifiers(const gvcf_options& opt,
                   const gvcf_deriv_options& dopt,
                   site_info& si,
                   calibration_models& model) {

    si.smod.clear();

    si.smod.is_unknown=(si.ref=='N');

    si.smod.is_used_covered=(si.n_used_calls!=0);
    si.smod.is_covered=(si.smod.is_used_covered || si.n_unused_calls!=0);

    if     (si.smod.is_unknown) {
        si.smod.gqx=0;
        si.smod.gq=0;
        si.smod.max_gt=0;
    } else if (si.dgt.genome.max_gt != si.dgt.poly.max_gt) {
        si.smod.gqx=0;
        si.smod.gq=si.dgt.poly.max_gt_qphred;
        si.smod.max_gt=si.dgt.poly.max_gt;
    } else {
        if (si.dgt.genome.max_gt_qphred<si.dgt.poly.max_gt_qphred) {
            set_site_gt(si.dgt.genome,si.smod);
        } else {
            set_site_gt(si.dgt.poly,si.smod);
        }
        si.smod.gq=si.dgt.poly.max_gt_qphred;
    }

    set_site_filters_CM(opt,dopt,si,model);
}



gvcf_aggregator::
gvcf_aggregator(const starling_options& opt,
                const starling_deriv_options& dopt,
                const reference_contig_segment& ref,
                std::ostream* osptr)
    : _opt(opt)
    , _report_range(dopt.report_range.begin_pos,dopt.report_range.end_pos)
    , _ref(ref)
    , _osptr(osptr)
    , _chrom(opt.bam_seq_name.c_str())
    , _indel_end_pos(0)
    , _indel_buffer_size(0)
    , _site_buffer_size(0)
    , _block(_opt.gvcf)
    , _head_pos(dopt.report_range.begin_pos)
{
    assert(_report_range.is_begin_pos);
    assert(_report_range.is_end_pos);

    if (! opt.is_gvcf_output()) return;

    // initialize codonPhaser
    if(_opt.do_codon_phasing){
        codon_phaser = Codon_phaser();
//        codon_phaser = Codon_phaser();
        #ifdef DEBUG_GVCF
            //log_os << "I have a phaser" << "\n";
        #endif
    }
    // initialize calibration model
    this->CM.load_models(opt.calibration_models_filename);
    this->CM.set_model(opt.calibration_model);


    assert(NULL != _osptr);
    assert((NULL !=_chrom) && (strlen(_chrom)>0));

    cdmap_t chrom_depth;
    if (_opt.gvcf.is_max_depth_factor && (! _opt.gvcf.chrom_depth_file.empty())) {
        parse_chrom_depth(_opt.gvcf.chrom_depth_file,chrom_depth);
        //TODO, verify that chroms match bam chroms

        cdmap_t::const_iterator cdi(chrom_depth.find(std::string(_chrom)));
        if (cdi == chrom_depth.end()) {
            std::ostringstream oss;
            oss << "ERROR: Can't find chromosome: '" << _chrom << "' in chrom depth file: " << _opt.gvcf.chrom_depth_file << "\n";
            throw blt_exception(oss.str().c_str());
        }
        _dopt.max_depth=(cdi->second*_opt.gvcf.max_depth_factor);
        assert(_dopt.max_depth>=0.);
        _dopt.is_max_depth=true;
    }

    if (! _opt.gvcf.is_skip_header) {
        finish_gvcf_header(_opt.gvcf,chrom_depth,dopt.bam_header_data,*_osptr);
    }

    add_site_modifiers(_opt.gvcf,_dopt,_empty_site,this->CM);
}



gvcf_aggregator::
~gvcf_aggregator() { flush(); }



// fill in missing sites:
void
gvcf_aggregator::
skip_to_pos(const pos_t target_pos) {

    // advance through any indel region by adding individual sites
    while (_head_pos<target_pos) {
        add_site_internal(get_empty_site(_head_pos));

        // only add one empty site after completing any pre-existing indel blocks,
        // then extend the block size of that one site as required:
        if (0 != _indel_buffer_size) continue;
        if (_opt.gvcf.is_block_compression) {
            assert(_block.count!=0);
            _block.count += (target_pos-_head_pos);
            _head_pos=target_pos;
        } else {
            _head_pos++;
        }
    }
}


void
gvcf_aggregator::
add_site(site_info& si) {

    add_site_modifiers(_opt.gvcf,_dopt,si,this->CM);
    skip_to_pos(si.pos);
    if (_opt.do_codon_phasing){
        bool emptyBuffer = codon_phaser.add_site(si);

        // Was site absorbed, if not release all buffered sites
        // and add these into the gVCF pipeline
        if (!codon_phaser.is_in_block || emptyBuffer){
//            codon_phaser.write_out_buffer();
//            log_os << "### buffer size:" << codon_phaser.buffer.size() << "\n";
            for (std::vector<site_info>::iterator it = codon_phaser.buffer.begin(); it != codon_phaser.buffer.end(); ++it){
//                log_os << *it << "\n";
                add_site_internal(*it);
            }
            skip_to_pos(si.pos);
            add_site_internal(si);
            codon_phaser.clear_buffer();
        }

    }
    else{
        add_site_internal(si);
    }
}

//Add sites to queue for writing to gVCF

void
gvcf_aggregator::
add_site_internal(const site_info& si) {

    _head_pos=si.pos+1;

    // resolve any current or previous indels before queue-ing site:
    if (0 != _indel_buffer_size) {
        if (si.pos>=_indel_end_pos) {
            process_overlaps();
        } else {
            while (_site_buffer.size() <= _site_buffer_size) {
                _site_buffer.push_back(site_info());
            }
            _site_buffer[_site_buffer_size++] = si;
            return;
        }
    }

    // write_site:
    queue_site_record(si);
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
    if (ik.is_breakpoint()) return;

    // don't handle homozygous reference calls for now:
    if (is_no_indel(dindel)) return;

    skip_to_pos(pos);

    // check if an indel is already buffered and we done't overlap it,
    // in which case we need to clear it first -- note this definition
    // of overlap deliberately picks up adjacent deletions:
    if ((0 != _indel_buffer_size) && (pos>_indel_end_pos)) {
        process_overlaps();
    }

    while (_indel_buffer.size() <= _indel_buffer_size) {
        _indel_buffer.push_back(indel_info());
    }
    _indel_buffer[_indel_buffer_size++].init(pos,ik,dindel,iri,isri);
    _indel_end_pos=std::max(_indel_end_pos,ik.right_pos());
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
get_hap_cigar(ALIGNPATH::path_t& apath,
              const indel_key& ik,
              const unsigned lead=1,
              const unsigned trail=0) {

    using namespace ALIGNPATH;

    apath.clear();
    if (lead) {
        apath.push_back(path_segment(MATCH,lead));
    }
    if (ik.delete_length()) {
        apath.push_back(path_segment(DELETE,ik.delete_length()));
    }
    if (ik.insert_length()) {
        apath.push_back(path_segment(INSERT,ik.insert_length()));
    }
    if (trail) {
        apath.push_back(path_segment(MATCH,trail));
    }
}


// figure out the per-site ploidy inside of indel based on each haplotype's match descriptor:
static
void
add_cigar_to_ploidy(const ALIGNPATH::path_t& apath,
                    std::vector<unsigned>& ploidy) {

    using namespace ALIGNPATH;
    int offset(-1);
    BOOST_FOREACH(const path_segment& ps, apath) {
        if (ps.type==MATCH) {
            for (unsigned j(0); j<ps.length; ++j) {
                if (offset>=0) ploidy[offset]++;
                offset++;
            }
        } else if (ps.type==DELETE) {
            offset+=ps.length;
        }
    }
}




static
void
add_indel_modifiers(const gvcf_options& opt,
                    const gvcf_deriv_options& dopt,
                    indel_info& ii) {

    if (ii.dindel.max_gt != ii.dindel.max_gt_poly) {
        ii.imod.gqx=0;
    } else {
        ii.imod.gqx=std::min(ii.dindel.max_gt_poly_qphred,ii.dindel.max_gt_qphred);
    }
    ii.imod.max_gt=ii.dindel.max_gt_poly;
    ii.imod.gq=ii.dindel.max_gt_poly_qphred;


    if (opt.is_min_gqx) {
        if (ii.imod.gqx<opt.min_gqx) ii.imod.set_filter(VCF_FILTERS::LowGQX);
    }

    if (dopt.is_max_depth) {
        if (ii.isri.depth > dopt.max_depth) ii.imod.set_filter(VCF_FILTERS::HighDepth);
    }

    if (opt.is_max_ref_rep) {
        if (ii.iri.is_repeat_unit) {
            if ((ii.iri.repeat_unit.size() <= 2) &&
                (static_cast<int>(ii.iri.ref_repeat_count) > opt.max_ref_rep)) {
                ii.imod.set_filter(VCF_FILTERS::HighRefRep);
            }
        }
    }
}





// is the current site eligible to even be considered for block compression?
static
bool
is_site_record_blockable(const gvcf_options& opt,
                         const site_info& si) {

    if (! opt.is_block_compression) return false;

    if (si.dgt.is_snp) return false;

    if (si.ref!='N') {
        const double reffrac(static_cast<double>(si.known_counts[si.dgt.ref_gt]) /
                             static_cast<double>(si.n_used_calls));
        if (reffrac+opt.block_max_nonref <= 1) return false;
    }
    return true;
}



// queue site record for writing, after
// possibly joining it into a compressed non-variant block
//
void
gvcf_aggregator::
queue_site_record(const site_info& si) {

    if (! is_site_record_blockable(_opt.gvcf,si)) {
        write_block_site_record();
        write_site_record(si);
        return;
    }

    if (! _block.test(si)) {
        write_block_site_record();
    }

    _block.join(si);
}



static
void
print_vcf_alt(const unsigned gt,
              const unsigned ref_gt,
              std::ostream& os) {

    bool is_print(false);
    for (unsigned b(0); b<N_BASE; ++b) {
        if (b==ref_gt) continue;
        if (DIGT::expect2(b,gt)) {
            if (is_print) os << ',';
            os << id_to_base(b);
            is_print=true;
        }
    }
    if (! is_print) os << '.';
}



static
void
print_site_ad(const site_info& si,
              std::ostream& os) {

    os << si.known_counts[si.dgt.ref_gt];

    for (unsigned b(0); b<N_BASE; ++b) {
        if (b==si.dgt.ref_gt) continue;
        if (DIGT::expect2(b,si.smod.max_gt)) {
            os << ',' << si.known_counts[b];
        }
    }
}


//writes out a SNP or block record
void
gvcf_aggregator::
write_site_record(const site_info& si) const {

    std::ostream& os(*_osptr);

    os << _chrom << '\t'  // CHROM
       << (si.pos+1) << '\t'  // POS
       << ".\t"           // ID
       << si.ref << '\t'; // REF

    // ALT
    if (si.smod.is_unknown || si.smod.is_block) {
        os << '.';
    } else {
        print_vcf_alt(si.smod.max_gt,si.dgt.ref_gt,os);
    }
    os << '\t';

    // QUAL:
    if (si.is_qual()) {
        os << si.dgt.genome.snp_qphred;
    } else {
        os << '.';
    }
    os << '\t';

    // FILTER:
    si.smod.write_filters(os);
    os << '\t';

    // INFO:
    if (si.smod.is_block) {
        if (_block.count>1) {
            os << "END=" << (si.pos+_block.count) << ';';
            os << _opt.gvcf.block_label;
        } else {
            os << '.';
        }
    } else {
        if (si.dgt.is_snp) {
            os << "SNVSB=";
            std::ofstream tmp_os;
            tmp_os.copyfmt(os);
            os << std::fixed << std::setprecision(1) << si.dgt.sb;
            os.copyfmt(tmp_os);
            os << ';';
            os << "SNVHPOL=" << si.hpol;
            if (_opt.is_compute_hapscore) {
                os << ';';
                os << "HaplotypeScore=" << si.hapscore;
            }
            //reported q-score
            if (si.Qscore>0){
                os << ';';
                os << "Qscore=" << si.Qscore;
            }

            // compute metrics nessacery
            if (_opt.is_compute_VQSRmetrics) {
                os << ';';
                os << "MQ=" << si.MQ;

                //if we have a het, report these metrics as well
//                if(si.get_gt()=="0/1"){
                os << ';';
                os << "MQRankSum=" << si.MQRankSum;
                os << ';';
                os << "BaseQRankSum=" << si.BaseQRankSum;
                os << ';';
                os << "ReadPosRankSum=" << si.ReadPosRankSum;
                os << ';';
                os << "DP=" << (si.n_used_calls+si.n_unused_calls);
                os << ';';
                os << "GQ=" << si.smod.gq;
                os << ';';
                os << "GQX=" << si.smod.gqx;
//                }
            }

        } else {
            os << '.';
        }
    }
    os << '\t';

    //FORMAT
    os << "GT";
    if (si.dgt.is_snp) {
        os << ":GQ";
    }
    os << ":GQX:DP:DPF";
    if (! si.smod.is_block) {
        os << ":AD";
    }
    os << '\t';

    //SAMPLE
    os << si.get_gt() << ':';
    if (si.dgt.is_snp) {
        os << si.smod.gq << ':';
    }
    if (si.smod.is_gqx()) {
        if (si.smod.is_block) {
            os << _block.block_gqx.min();
        } else {
            os << si.smod.gqx;
        }
    } else {
        os << '.';
    }
    os << ':';
    //print DP:DPF
    if (si.smod.is_block) {
        os << _block.block_dpu.min() << ':'
           << _block.block_dpf.min();
    } else {
        os << si.n_used_calls << ':'
           << si.n_unused_calls;
    }

    if (! si.smod.is_block) {
        // print AD
        os << ':';
        print_site_ad(si,os);
    }

    os << '\n';
}


// set the CIGAR string:
void
gvcf_aggregator::
modify_single_indel_record() {
    assert(_indel_buffer_size==1);

    indel_info& ii(_indel_buffer[0]);
    get_hap_cigar(ii.imod.cigar,ii.ik);

    add_indel_modifiers(_opt.gvcf,_dopt,ii);
}



static
void
modify_indel_overlap_site(const gvcf_options& opt,
                          const gvcf_deriv_options& dopt,
                          const indel_info& ii,
                          const unsigned ploidy,
                          site_info& si) {

#ifdef DEBUG_GVCF
    log_os << "CHIRP: indel_overlap_site smod before: " << si.smod << "\n";
    log_os << "CHIRP: indel_overlap_site imod before: " << ii.imod << "\n";
#endif

    // inherit any filters from the indel:
    si.smod.filters |= ii.imod.filters;

#ifdef DEBUG_GVCF
    log_os << "CHIRP: indel_overlap_site smod after: " << si.smod << "\n";
#endif

    // limit qual and gq values to those of the indel
    si.dgt.genome.snp_qphred = std::min(si.dgt.genome.snp_qphred,ii.dindel.indel_qphred);
    si.smod.gqx = std::min(si.smod.gqx,ii.dindel.max_gt_qphred);

    // change ploidy:
    if (ploidy==1) {
        if (DIGT::is_het(si.smod.max_gt)) {
            si.smod.set_filter(VCF_FILTERS::SiteConflict);
            //si.smod.modified_gt=MODIFIED_SITE_GT::UNKNOWN;
        } else {
            if (si.smod.max_gt == si.dgt.ref_gt) {
                si.smod.modified_gt=MODIFIED_SITE_GT::ZERO;
            } else {
                si.smod.modified_gt=MODIFIED_SITE_GT::ONE;
            }
        }
    } else if (ploidy==0) {
        if (si.smod.max_gt == si.dgt.ref_gt) {
            si.smod.modified_gt=MODIFIED_SITE_GT::UNKNOWN;
            si.smod.is_zero_ploidy=true;
        } else {
            si.smod.set_filter(VCF_FILTERS::SiteConflict);
        }
    } else {
        assert(0);
    }

    // after all those changes we need to rerun the site filters:
    set_site_filters(opt,dopt,si);
}



static
void
modify_indel_conflict_site(site_info& si) {

    si.smod.set_filter(VCF_FILTERS::IndelConflict);
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

    ii.imod.ploidy.resize(_indel_end_pos-ii.pos,0);

    // add per-haplotype information:
    for (unsigned hap(0); hap<2; ++hap) {
        //reduce qual and gt to the lowest of the set:
        if (hap) {
            if (ii.dindel.indel_qphred>_indel_buffer[hap].dindel.indel_qphred) {
                ii.dindel.indel_qphred = _indel_buffer[hap].dindel.indel_qphred;
            }
            if (ii.dindel.max_gt_qphred>_indel_buffer[hap].dindel.max_gt_qphred) {
                ii.dindel.max_gt_qphred = _indel_buffer[hap].dindel.max_gt_qphred;
            }
        }

        // extend leading sequence start back 1 for vcf compat, and end back 1 to concat with vcf_indel_seq
        _ref.get_substring(indel_begin_pos,(_indel_buffer[hap].pos-indel_begin_pos)-1,leading_seq);
        const unsigned trail_len(_indel_end_pos-_indel_buffer[hap].ik.right_pos());
        _ref.get_substring(_indel_end_pos-trail_len,trail_len,trailing_seq);


        _indel_buffer[hap].iri.vcf_indel_seq = leading_seq + _indel_buffer[hap].iri.vcf_indel_seq + trailing_seq;

        get_hap_cigar(_indel_buffer[hap].imod.cigar,
                      _indel_buffer[hap].ik,
                      leading_seq.size()+1,
                      trailing_seq.size());

        // add to the ploidy object:
        add_cigar_to_ploidy(_indel_buffer[hap].imod.cigar,ii.imod.ploidy);

        add_indel_modifiers(_opt.gvcf,_dopt,_indel_buffer[hap]);
        if (hap>0) {
            ii.imod.filters |= _indel_buffer[hap].imod.filters;
        }
    }
}



// set the CIGAR string:
void
gvcf_aggregator::
modify_conflict_indel_record() {
    assert(_indel_buffer_size>1);

    for (unsigned i(0); i<_indel_buffer_size; ++i) {
        indel_info& ii(_indel_buffer[i]);
        get_hap_cigar(ii.imod.cigar,ii.ik);

        ii.imod.set_filter(VCF_FILTERS::IndelConflict);

        add_indel_modifiers(_opt.gvcf,_dopt,ii);
    }
}



void
gvcf_aggregator::
write_indel_record(const unsigned write_index) {

    assert(_indel_buffer_size>0);

    // flush any non-variant block before starting:
    write_block_site_record();

    std::ostream& os(*_osptr);
    indel_info& ii(_indel_buffer[write_index]);

    os << _chrom << '\t'   // CHROM
       << ii.pos << '\t'   // POS
       << ".\t"            // ID
       << ii.iri.vcf_ref_seq << '\t'; // REF

    // ALT
    unsigned end_index(write_index);
    if (ii.imod.is_overlap) {
        end_index++;
    }

    for (unsigned i(write_index); i<=end_index; ++i) {
        if (i!=write_index) os << ',';
        os << _indel_buffer[i].iri.vcf_indel_seq;
    }
    os << '\t';

    os << ii.dindel.indel_qphred << '\t'; //QUAL

    // FILTER:
    ii.imod.write_filters(os);
    os << '\t';

    // INFO
    os << "CIGAR=";
    for (unsigned i(write_index); i<=end_index; ++i) {
        if (i!=write_index) os << ',';
        os << _indel_buffer[i].imod.cigar;
    }
    os << ';';
    os << "RU=";
    for (unsigned i(write_index); i<=end_index; ++i) {
        if (i!=write_index) os << ',';
        if (_indel_buffer[i].iri.is_repeat_unit &&
            (_indel_buffer[i].iri.repeat_unit.size() <= 20)) {
            os << _indel_buffer[i].iri.repeat_unit;
        } else {
            os << '.';
        }
    }
    os << ';';
    os << "REFREP=";
    for (unsigned i(write_index); i<=end_index; ++i) {
        if (i!=write_index) os << ',';
        if (_indel_buffer[i].iri.is_repeat_unit) {
            os << _indel_buffer[i].iri.ref_repeat_count;
        } else {
            os << '.';
        }
    }
    os << ';';
    os << "IDREP=";
    for (unsigned i(write_index); i<=end_index; ++i) {
        if (i!=write_index) os << ',';
        if (_indel_buffer[i].iri.is_repeat_unit) {
            os << _indel_buffer[i].iri.indel_repeat_count;
        } else {
            os << '.';
        }
    }

    if (_opt.is_compute_VQSRmetrics)
    {
        os << ';';
        os << "MQ=" << ii.MQ;

        //if we have a het, report these metrics as well
        //                if(si.get_gt()=="0/1"){
        os << ';';
        os << "MQRankSum=" << ii.MQRankSum;
        os << ';';
        os << "BaseQRankSum=" << ii.BaseQRankSum;
        os << ';';
        os << "ReadPosRankSum=" << ii.ReadPosRankSum;
    }

    //write out Q-score
    //if not ii.

    os << '\t';


    //FORMAT
    os << "GT:GQ:GQX:DPI:AD" << '\t';

    //SAMPLE
    os << ii.get_gt() << ':'
       << ii.imod.gq << ':'
       << ii.imod.gqx  << ':'
       << ii.isri.depth << ':';

    // SAMPLE AD:
    unsigned ref_count(0);
    for (unsigned i(write_index); i<=end_index; ++i) {
        ref_count = std::max(ref_count,_indel_buffer[i].isri.n_q30_ref_reads);
    }
    os << ref_count;
    for (unsigned i(write_index); i<=end_index; ++i) {
        os << ',' << _indel_buffer[i].isri.n_q30_indel_reads;
    }

    os << '\n';
}



void
gvcf_aggregator::
process_overlaps() {

    if (0==_indel_buffer_size) return;

    bool is_conflict_print(false);

    // do the overlap processing:
    if (_indel_buffer_size==1) {
        // simple case of no overlap:
        modify_single_indel_record();
    } else {
        if (is_simple_indel_overlap(_indel_buffer,_indel_buffer_size)) {
            // handle the simplest possible overlap case (two hets):
            modify_overlap_indel_record();
        } else {
            // mark the whole region as conflicting
            modify_conflict_indel_record();
            is_conflict_print=true;
        }
    }

    //    *_osptr << "INDEL_SIZE: " << _indel_buffer_size << "\n";

    // process sites to be consistent with overlapping indels:
    for (unsigned i(0); i<_site_buffer_size; ++i) {

#ifdef DEBUG_GVCF
        log_os << "CHIRP: indel overlapping site: " << _site_buffer[i].pos << "\n";
#endif
        const pos_t offset(_site_buffer[i].pos-_indel_buffer[0].pos);
        assert(offset>=0);
        if (! is_conflict_print) {
            modify_indel_overlap_site(_opt.gvcf,_dopt,
                                      _indel_buffer[0],
                                      _indel_buffer[0].get_ploidy(offset),
                                      _site_buffer[i]);
        } else {
            modify_indel_conflict_site(_site_buffer[i]);
        }
    }


    unsigned indel_index(0);
    unsigned site_index(0);

    while (true) {
        const bool is_indel(indel_index<_indel_buffer_size);
        const bool is_site(site_index<_site_buffer_size);
        if (! (is_indel || is_site)) break;

        if (is_indel && ((! is_site) || _indel_buffer[indel_index].pos <= _site_buffer[site_index].pos)) {
            // print indel:
            write_indel_record(indel_index);
            if (is_conflict_print) {
                indel_index++;
            } else {
                indel_index=_indel_buffer_size;
            }
        } else {
            // print site:
            log_os << "site record" << "\n";
            queue_site_record(_site_buffer[site_index]);
            site_index++;
        }
    }
    _indel_buffer_size = 0;
    _site_buffer_size = 0;
}




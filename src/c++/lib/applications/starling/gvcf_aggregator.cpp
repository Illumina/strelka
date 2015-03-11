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

#include "gvcf_aggregator.hh"
#include "gvcf_header.hh"
#include "blt_util/blt_exception.hh"
#include "blt_util/chrom_depth_map.hh"
#include "blt_util/io_util.hh"

#include <iomanip>
#include <iostream>
#include <sstream>


//#define DEBUG_GVCF


#ifdef DEBUG_GVCF
#include "blt_util/log.hh"
#endif



static
void
set_site_gt(const diploid_genotype::result_set& rs,
            site_modifiers& smod)
{
    smod.max_gt=rs.max_gt;
    smod.gqx=rs.max_gt_qphred;
    smod.gq  = 2;
}



static
void
add_site_modifiers(site_info& si,
                   calibration_models& model)
{
    si.smod.clear();
    si.smod.is_unknown=(si.ref=='N');
    si.smod.is_used_covered=(si.n_used_calls!=0);
    si.smod.is_covered=(si.smod.is_used_covered || si.n_unused_calls!=0);

    if     (si.smod.is_unknown)
    {
        si.smod.gqx=0;
        si.smod.gq=0;
        si.smod.max_gt=0;
    }
    else if (si.dgt.genome.max_gt != si.dgt.poly.max_gt)
    {
        si.smod.gqx=0;
        si.smod.gq=si.dgt.poly.max_gt_qphred;
        si.smod.max_gt=si.dgt.poly.max_gt;
    }
    else
    {
        if (si.dgt.genome.max_gt_qphred<si.dgt.poly.max_gt_qphred)
        {
            set_site_gt(si.dgt.genome,si.smod);
        }
        else
        {
            set_site_gt(si.dgt.poly,si.smod);
        }
        si.smod.gq=si.dgt.poly.max_gt_qphred;
    }

    model.clasify_site(si);
}



void gvcf_aggregator::write_block_site_record()
{
    if (_block.count<=0) return;
    write_site_record(_block.record);
    _block.reset();
}

gvcf_aggregator::
gvcf_aggregator(
    const starling_options& opt,
    const starling_deriv_options& dopt,
    const reference_contig_segment& ref,
    const RegionTracker& nocompress_regions,
    std::ostream* osptr,
    starling_read_buffer& read_buffer,
    const unsigned max_read_len)
    : _opt(opt)
    , _report_range(dopt.report_range.begin_pos,dopt.report_range.end_pos)
    , _ref(ref)
    , _osptr(osptr)
    , _chrom(opt.bam_seq_name.c_str())
    , _dopt(dopt.gvcf)
    , _indel_end_pos(0)
    , _indel_buffer_size(0)
    , _site_buffer_size(0)
    , _block(_opt.gvcf)
    , _head_pos(dopt.report_range.begin_pos)
    , _CM(_opt, dopt.gvcf)
    , _gvcf_comp(opt.gvcf,nocompress_regions)
    , _codon_phaser(opt, read_buffer, max_read_len)
//    , _assemble_stream(opt, read_buffer, max_read_len,nocompress_regions)
{
    assert(_report_range.is_begin_pos);
    assert(_report_range.is_end_pos);

    if (! opt.gvcf.is_gvcf_output()) return;

    assert(nullptr != _osptr);
    assert((nullptr !=_chrom) && (strlen(_chrom)>0));

    if (! _opt.gvcf.is_skip_header)
    {
        finish_gvcf_header(_opt,_dopt, _dopt.chrom_depth,dopt.bam_header_data,*_osptr,this->_CM);
    }

    add_site_modifiers(_empty_site,this->_CM);
}

gvcf_aggregator::
~gvcf_aggregator()
{
    flush();
}

bool
gvcf_aggregator::
add_site(site_info& si)
{
    add_site_modifiers(si, _CM);

    if (si.dgt.is_haploid())
    {
        if (si.smod.max_gt == si.dgt.ref_gt)
        {
            si.smod.modified_gt=MODIFIED_SITE_GT::ZERO;
        }
        else
        {
            si.smod.modified_gt=MODIFIED_SITE_GT::ONE;
        }
    }
    else if (si.dgt.is_noploid())
    {
        si.smod.set_filter(VCF_FILTERS::PloidyConflict);
    }

//    if (_opt.do_codon_phasing
//        && (si.is_het() || _codon_phaser.is_in_block()))
//    {
//        const bool emptyBuffer = _codon_phaser.add_site(si);
//        if (!_codon_phaser.is_in_block() || emptyBuffer)
//            this->output_phased_blocked();
//    }
//    else if (_opt.do_assemble)
//    {
//    	const bool emptyBuffer = _assembler.add_site(si);
//        if (emptyBuffer)
//        	this->output_phased_blocked();
//    }

    else
    {
        skip_to_pos(si.pos);
        add_site_internal(si);
    }
    return true;
}


// fill in missing sites
void
gvcf_aggregator::
skip_to_pos(const pos_t target_pos)
{
    // advance through any indel region by adding individual sites
    while (_head_pos<target_pos)
    {
        const site_info& si = get_empty_site(_head_pos);
        add_site_internal(si);
        // only add one empty site after completing any pre-existing indel blocks,
        // then extend the block size of that one site as required:
        if (0 != _indel_buffer_size) continue;

        if (_gvcf_comp.is_range_compressable(known_pos_range2(si.pos,target_pos)))
        {
            assert(_block.count!=0);
            _block.count += (target_pos-_head_pos);
            _head_pos= target_pos;
        }
    }
}



void
gvcf_aggregator::
output_phased_blocked()
{
    // output the codon-phaser or assembler buffer to gVCF queue

    //case assembler
//	if (_opt.do_assemble){
//		for (const site_info& si : _assembler.buffer())
//			{
//				this->skip_to_pos(si.pos);
//				add_site_internal(si);
//			}
//			_assembler.clear();
//	}
    // case codon-phaser
//    for (const site_info& si : _codon_phaser.buffer())
//    {
//        this->skip_to_pos(si.pos);
//        add_site_internal(si);
//    }
//    _codon_phaser.clear();
}

//Add sites to queue for writing to gVCF
void
gvcf_aggregator::
add_site_internal(const site_info& si)
{
    if (si.smod.is_phased_region)
        _head_pos=si.pos+si.phased_ref.length();
    else
        _head_pos=si.pos+1;

    // resolve any current or previous indels before queue-ing site:
    if (0 != _indel_buffer_size)
    {
        if (si.pos>=_indel_end_pos)
        {
            process_overlaps();
        }
        else
        {
            while (_site_buffer.size() <= _site_buffer_size)
            {
                _site_buffer.emplace_back();
            }
            _site_buffer[_site_buffer_size++] = si;
            return;
        }
    }

    // write_site
    queue_site_record(si);
}

static
bool
is_het_indel(const starling_diploid_indel_core& dindel)
{
    return (dindel.max_gt==STAR_DIINDEL::HET);
}

static
bool
is_no_indel(const starling_diploid_indel_core& dindel)
{
    return (dindel.max_gt==STAR_DIINDEL::NOINDEL);
}

bool
gvcf_aggregator::
add_indel(const pos_t pos,
          const indel_key ik,
          const starling_diploid_indel_core& dindel,
          const starling_indel_report_info& iri,
          const starling_indel_sample_report_info& isri)
{
    // we can't handle breakends at all right now:
    if (ik.is_breakpoint()) return true;

    // don't handle homozygous reference calls unless genotyping is forced
    if (is_no_indel(dindel) && !dindel.is_forced_output) return true;

    // if we are in phasing a block and encounter an indel, make sure we empty block before doing anything else
    if (_opt.do_codon_phasing && this->_codon_phaser.is_in_block())
        this->output_phased_blocked();

    skip_to_pos(pos);

    // check if an indel is already buffered and
    // either we don't overlap it or we get homRef for the forced-genotyped indel,
    // in which case we need to clear it first -- note this definition
    // of overlap deliberately picks up adjacent deletions:
    if ((0 != _indel_buffer_size) && ((pos>_indel_end_pos) || is_no_indel(dindel)))
    {
        process_overlaps();
    }

    while (_indel_buffer.size() <= _indel_buffer_size)
    {
        _indel_buffer.emplace_back();
    }
    _indel_buffer[_indel_buffer_size++].init(pos,ik,dindel,iri,isri);
    _indel_end_pos=std::max(_indel_end_pos,ik.right_pos());

    // add filter for all indels in no-ploid regions:
    if (dindel.is_noploid())
    {
        indel_info& ii(_indel_buffer[_indel_buffer_size-1]);
        ii.imod.set_filter(VCF_FILTERS::PloidyConflict);
    }

    // clear the current homRef indel
    if (is_no_indel(dindel))
    {
        process_overlaps();
    }
    return true;
}

static
bool
is_simple_indel_overlap(const std::vector<indel_info>& indel_buffer,
                        const unsigned size)
{
    return (size==2 &&
            is_het_indel(indel_buffer[0].dindel) &&
            is_het_indel(indel_buffer[1].dindel));
}


static
void
get_hap_cigar(ALIGNPATH::path_t& apath,
              const indel_key& ik,
              const unsigned lead=1,
              const unsigned trail=0)
{
    using namespace ALIGNPATH;

    apath.clear();
    if (lead)
    {
        apath.push_back(path_segment(MATCH,lead));
    }
    if (ik.delete_length())
    {
        apath.push_back(path_segment(DELETE,ik.delete_length()));
    }
    if (ik.insert_length())
    {
        apath.push_back(path_segment(INSERT,ik.insert_length()));
    }
    if (trail)
    {
        apath.push_back(path_segment(MATCH,trail));
    }
}

// figure out the per-site ploidy inside of indel based on each haplotype's match descriptor:
static
void
add_cigar_to_ploidy(const ALIGNPATH::path_t& apath,
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

// queue site record for writing, after
// possibly joining it into a compressed non-variant block
//
void
gvcf_aggregator::
queue_site_record(const site_info& si)
{
    //test for basic blocking criteria
    if (! _gvcf_comp.is_site_compressable(si))
    {
        write_block_site_record();
        write_site_record(si);
        return;
    }

    if (! _block.test(si))
    {
        write_block_site_record();
    }
    _block.join(si);
}



static
void
print_vcf_alt(const unsigned gt,
              const unsigned ref_gt,
              std::ostream& os)
{
    bool is_print(false);
    for (unsigned b(0); b<N_BASE; ++b)
    {
        if (b==ref_gt) continue;
        if (DIGT::expect2(b,gt))
        {
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
              std::ostream& os)
{
    os << si.known_counts[si.dgt.ref_gt];

    for (unsigned b(0); b<N_BASE; ++b)
    {
        if (b==si.dgt.ref_gt) continue;
        if (DIGT::expect2(b,si.smod.max_gt))
        {
            os << ',' << si.known_counts[b];
        }
    }
}


//writes out a SNP or block record
void
gvcf_aggregator::
write_site_record(const site_info& si) const
{
    std::ostream& os(*_osptr);

    os << _chrom << '\t'  // CHROM
       << (si.pos+1) << '\t'  // POS
       << ".\t";           // ID

    if (si.smod.is_phased_region)
        os  << si.phased_ref << '\t'; // REF
    else
        os  << si.ref << '\t'; // REF

    // ALT
    if ((si.smod.is_unknown || si.smod.is_block) && !si.smod.is_phased_region)
    {
        os << '.';
    }
    else
    {
        if (si.smod.is_phased_region)
        {
            os << si.phased_alt;

        }
        else
            print_vcf_alt(si.smod.max_gt,si.dgt.ref_gt,os);
    }
    os << '\t';

    // QUAL:
    if (si.is_qual())
    {
        os << si.dgt.genome.snp_qphred;
    }
    else
    {
        os << '.';
    }
    os << '\t';

    // FILTER:
    si.smod.write_filters(os);
    os << '\t';

    // INFO:
    if (si.smod.is_block)
    {
        if (_block.count>1)
        {
            os << "END=" << (si.pos+_block.count) << ';';
            os << _dopt.block_label;
        }
        else
        {
            os << '.';
        }
    }
    else
    {
        if (si.dgt.is_snp)
        {
            os << "SNVSB=";
            {
                const StreamScoper ss(os);
                os << std::fixed << std::setprecision(1) << si.dgt.sb;
            }
            os << ';';
            os << "SNVHPOL=" << si.hpol;
            if (_opt.is_compute_hapscore)
            {
                os << ';';
                os << "HaplotypeScore=" << si.hapscore;
            }

            // TODO only report VQSR metrics if flag explicitly set, not for calibration model
            if (_opt.is_report_germline_VQSRmetrics)
            {
                os << ';';
                os << "MQ=" << si.MQ;

                os << ';';
                os << "MQRankSum=" << si.MQRankSum;
                os << ';';
                os << "BaseQRankSum=" << si.BaseQRankSum;
                os << ';';
                os << "ReadPosRankSum=" << si.ReadPosRankSum;
                os << ';';
                os << "AvgBaseQ=" << si.avgBaseQ;
                os << ';';
                os << "AvgPos=" << si.rawPos;
// if you uncomment the following, make sure you also uncomment the matching INFO header entry in gvcf_header.cpp
//                os << ';';
//                os << "MapQ0Count=" << si.mapq_zero;

// N.B. DP is in FORMAT already, and that seems to be where Nondas's code expects to find it, so suppress it here:
//                os << ';';
//                os << "DP=" << (si.n_used_calls+si.n_unused_calls);

            }
//            //reported q-score
//            if (si.Qscore>0) {
//                os << ';';
//                os << "Qscore=" << si.Qscore;
//            }

        }
        else
        {
            os << '.';
        }
    }
    os << '\t';

    //FORMAT
    os << "GT";
    if (si.dgt.is_snp)
    {
        os << ":GQ";
    }
    os << ":GQX:DP:DPF";
    if (! si.smod.is_block || si.smod.is_phased_region)
    {
        os << ":AD";
    }
    os << '\t';

    //SAMPLE
    os << si.get_gt() << ':';
    if (si.dgt.is_snp)
    {
        os << si.smod.gq << ':';
    }
    if (si.smod.is_gqx())
    {
        if (si.smod.is_block && !si.smod.is_phased_region)
        {
            os << _block.block_gqx.min();
        }
        else
        {
            if (si.Qscore>=0)
                os << si.Qscore ;
            else
                os << si.smod.gqx;
        }
    }
    else
    {
        os << '.';
    }
    os << ':';
    //print DP:DPF
    if (si.smod.is_block && !si.smod.is_phased_region)
    {
        os << _block.block_dpu.min() << ':'
           << _block.block_dpf.min();
    }
    else
    {
        os << si.n_used_calls << ':'
           << si.n_unused_calls;
    }

    if (si.smod.is_phased_region)
        os << ':' << si.phased_AD;
    else if (! si.smod.is_block)
    {
        os << ':';
        print_site_ad(si,os);
    }
    os << '\n';
}

// set the CIGAR string:
void
gvcf_aggregator::
modify_single_indel_record()
{
    assert(_indel_buffer_size==1);

    indel_info& ii(_indel_buffer[0]);
    get_hap_cigar(ii.imod.cigar,ii.ik);

    _CM.clasify_site(ii);
}

static
void
modify_indel_overlap_site(const indel_info& ii,
                          const unsigned ploidy,
                          site_info& si,calibration_models& CM)
{
//#ifdef DEBUG_GVCF
//    log_os << "CHIRP: indel_overlap_site smod before: " << si.smod << "\n";
//    log_os << "CHIRP: indel_overlap_site imod before: " << ii.imod << "\n";
//#endif

    // inherit any filters from the indel:
    si.smod.filters |= ii.imod.filters;

#ifdef DEBUG_GVCF
    log_os << "CHIRP: indel_overlap_site smod after: " << si.smod << "\n";
#endif

    // limit qual and gq values to those of the indel
    si.dgt.genome.snp_qphred = std::min(si.dgt.genome.snp_qphred,ii.dindel.indel_qphred);
    si.smod.gqx = std::min(si.smod.gqx,ii.dindel.max_gt_qphred);

    // change ploidy:
    if (ploidy==1)
    {
        if (DIGT::is_het(si.smod.max_gt))
        {
            si.smod.set_filter(VCF_FILTERS::SiteConflict);
            //si.smod.modified_gt=MODIFIED_SITE_GT::UNKNOWN;
        }
        else
        {
            if (si.smod.max_gt == si.dgt.ref_gt)
            {
                si.smod.modified_gt=MODIFIED_SITE_GT::ZERO;
            }
            else
            {
                si.smod.modified_gt=MODIFIED_SITE_GT::ONE;
            }
        }
    }
    else if (ploidy==0)
    {
        if (si.smod.max_gt == si.dgt.ref_gt)
        {
            si.smod.modified_gt=MODIFIED_SITE_GT::UNKNOWN;
            si.smod.is_zero_ploidy=true;
        }
        else
        {
            si.smod.set_filter(VCF_FILTERS::SiteConflict);
        }
    }
    else if (ploidy!=2)
    {
        assert(0);
    }

    // after all those changes we need to rerun the site filters:
    CM.clasify_site(si);
}

static
void
modify_indel_conflict_site(site_info& si)
{
    si.smod.set_filter(VCF_FILTERS::IndelConflict);
}

void
gvcf_aggregator::
modify_overlap_indel_record()
{
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
    for (unsigned hap(0); hap<2; ++hap)
    {
        //reduce qual and gt to the lowest of the set:
        if (hap)
        {
            if (ii.dindel.indel_qphred>_indel_buffer[hap].dindel.indel_qphred)
            {
                ii.dindel.indel_qphred = _indel_buffer[hap].dindel.indel_qphred;
            }
            if (ii.dindel.max_gt_qphred>_indel_buffer[hap].dindel.max_gt_qphred)
            {
                ii.dindel.max_gt_qphred = _indel_buffer[hap].dindel.max_gt_qphred;
            }

            _indel_buffer[hap].imod.is_overlap=true;

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
        _CM.clasify_site(_indel_buffer[hap]);
        if (hap>0)
        {
            ii.imod.filters |= _indel_buffer[hap].imod.filters;
        }
    }
}



// set the CIGAR string:
void
gvcf_aggregator::
modify_conflict_indel_record()
{
    assert(_indel_buffer_size>1);

    for (unsigned i(0); i<_indel_buffer_size; ++i)
    {
        indel_info& ii(_indel_buffer[i]);
        get_hap_cigar(ii.imod.cigar,ii.ik);

        ii.imod.set_filter(VCF_FILTERS::IndelConflict);

        _CM.clasify_site(ii);
    }
}



void
gvcf_aggregator::
write_indel_record(const unsigned write_index)
{
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
    if (ii.imod.is_overlap)
    {
        end_index++;
    }

    for (unsigned i(write_index); i<=end_index; ++i)
    {
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
    for (unsigned i(write_index); i<=end_index; ++i)
    {
        if (i!=write_index) os << ',';
        os << _indel_buffer[i].imod.cigar;
    }
    os << ';';
    os << "RU=";
    for (unsigned i(write_index); i<=end_index; ++i)
    {
        if (i!=write_index) os << ',';
        if (_indel_buffer[i].iri.is_repeat_unit() &&
            (_indel_buffer[i].iri.repeat_unit.size() <= 20))
        {
            os << _indel_buffer[i].iri.repeat_unit;
        }
        else
        {
            os << '.';
        }
    }
    os << ';';
    os << "REFREP=";
    for (unsigned i(write_index); i<=end_index; ++i)
    {
        if (i!=write_index) os << ',';
        if (_indel_buffer[i].iri.is_repeat_unit())
        {
            os << _indel_buffer[i].iri.ref_repeat_count;
        }
        else
        {
            os << '.';
        }
    }
    os << ';';
    os << "IDREP=";
    for (unsigned i(write_index); i<=end_index; ++i)
    {
        if (i!=write_index) os << ',';
        if (_indel_buffer[i].iri.is_repeat_unit())
        {
            os << _indel_buffer[i].iri.indel_repeat_count;
        }
        else
        {
            os << '.';
        }
    }

//    if (ii.Qscore>0) {
//        os << ';';
//        os << "Qscore=" << ii.Qscore;
//    }

//    only report metrics if flag explicitly set
//    if (_opt.is_compute_VQSRmetrics)
//    {
//        os << ';';
//        os << "MQ=" << ii.MQ;
//        os << ';';
//        os << "MQRankSum=" << ii.MQRankSum;
//        os << ';';
//        os << "BaseQRankSum=" << ii.BaseQRankSum;
//        os << ';';
//        os << "ReadPosRankSum=" << ii.ReadPosRankSum;
//    }
    os << '\t';

    //FORMAT
    os << "GT:GQ:GQX:DPI:AD" << '\t';

    //SAMPLE
    os << ii.get_gt() << ':'
       << ii.imod.gq << ':';

    if (ii.Qscore>=0)
        os << ii.Qscore  << ':';
    else
        os << ii.imod.gqx  << ':';

    os << ii.isri.depth << ':';

    // SAMPLE AD:
    unsigned ref_count(0);
    for (unsigned i(write_index); i<=end_index; ++i)
    {
        ref_count = std::max(ref_count,_indel_buffer[i].isri.n_q30_ref_reads);
    }
    os << ref_count;
    for (unsigned i(write_index); i<=end_index; ++i)
    {
        os << ',' << _indel_buffer[i].isri.n_q30_indel_reads;
    }
    os << '\n';
}



void
gvcf_aggregator::
process_overlaps()
{
    if (0==_indel_buffer_size) return;

    bool is_conflict_print(false);

    // do the overlap processing:
    if (_indel_buffer_size==1)
    {
        // simple case of no overlap:
        modify_single_indel_record();
    }
    else
    {
        if (is_simple_indel_overlap(_indel_buffer,_indel_buffer_size))
        {
            // handle the simplest possible overlap case (two hets):
            modify_overlap_indel_record();
        }
        else
        {
            // mark the whole region as conflicting
            modify_conflict_indel_record();
            is_conflict_print=true;
        }
    }

    //    *_osptr << "INDEL_SIZE: " << _indel_buffer_size << "\n";

    // process sites to be consistent with overlapping indels:
    for (unsigned i(0); i<_site_buffer_size; ++i)
    {
#ifdef DEBUG_GVCF
        log_os << "CHIRP: indel overlapping site: " << _site_buffer[i].pos << "\n";
#endif
        const pos_t offset(_site_buffer[i].pos-_indel_buffer[0].pos);
        assert(offset>=0);
        if (! is_conflict_print)
        {
            modify_indel_overlap_site( _indel_buffer[0],
                                       _indel_buffer[0].get_ploidy(offset),
                                       _site_buffer[i], this->_CM);
        }
        else
        {
            modify_indel_conflict_site(_site_buffer[i]);
        }
    }

    unsigned indel_index(0);
    unsigned site_index(0);

    while (true)
    {
        const bool is_indel(indel_index<_indel_buffer_size);
        const bool is_site(site_index<_site_buffer_size);
        if (! (is_indel || is_site)) break;

        if (is_indel && ((! is_site) || _indel_buffer[indel_index].pos <= _site_buffer[site_index].pos))
        {
            // print indel:
            write_indel_record(indel_index);
            if (is_conflict_print)
            {
                indel_index++;
            }
            else
            {
                indel_index=_indel_buffer_size;
            }
        }
        else
        {
            // print site:
            //log_os << "site record" << "\n";
            queue_site_record(_site_buffer[site_index]);
            site_index++;
        }
    }
    _indel_buffer_size = 0;
    _site_buffer_size = 0;
}

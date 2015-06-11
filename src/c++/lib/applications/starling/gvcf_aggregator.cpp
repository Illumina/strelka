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
#include <stdexcept>


//#define DEBUG_GVCF


#ifdef DEBUG_GVCF
#include "blt_util/log.hh"
#endif




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
    const pos_basecall_buffer& bc_buff)
    : _opt(opt)
    , _report_range(dopt.report_range.begin_pos,dopt.report_range.end_pos)
    , _ref(ref)
    , _osptr(osptr)
    , _chrom(opt.bam_seq_name.c_str())
    , _dopt(dopt.gvcf)
    , _block(_opt.gvcf)
    , _head_pos(dopt.report_range.begin_pos)
    , _CM(_opt, dopt.gvcf)
    , _gvcf_comp(opt.gvcf,nocompress_regions)
    , _overlapper(_CM, _ref, *this)
    , _codon_phaser(opt, bc_buff, ref, _overlapper)
    , _targeted_region_processor(opt.gvcf.targeted_regions_bedfile, _chrom, _codon_phaser)
    , _head(_CM, _targeted_region_processor)
{
    assert(_report_range.is_begin_pos);
    assert(_report_range.is_end_pos);

    if (! opt.gvcf.is_gvcf_output())
        throw std::invalid_argument("gvcf_aggregator cannot be constructed with nothing to do.");

    assert(nullptr != _osptr);
    assert((nullptr !=_chrom) && (strlen(_chrom)>0));

    if (! _opt.gvcf.is_skip_header)
    {
        finish_gvcf_header(_opt,_dopt, _dopt.chrom_depth,dopt.bam_header_data, *_osptr, _CM);
    }

    variant_prefilter_stage::add_site_modifiers(_empty_site, _empty_site.smod, _CM);

}



gvcf_aggregator::
~gvcf_aggregator()
{
    _head.flush();
}

void
gvcf_aggregator::
add_site(
    site_info& si)
{
    _head.process(si);
}
void gvcf_aggregator::filter_site_by_last_indel_overlap(site_info& si)
{
    if (_last_indel)
    {
        if (si.pos >= _last_indel->end())
        {
            _last_indel.reset(nullptr);
        }
        else
        {
            _overlapper.modify_overlapping_site(*_last_indel, si);
        }
    }
}

// fill in missing sites
void
gvcf_aggregator::
skip_to_pos(const pos_t target_pos)
{
    // advance through any indel region by adding individual sites
    while (_head_pos<target_pos)
    {
        // TODO: this filtering should be delegated to the block compressor
        site_info si = get_empty_site(_head_pos);

        add_site_internal(si);
        // Don't do compressed ranges if there is an overlapping indel
		// filters are being applied to the overlapping positions
        if (_last_indel) continue;

        if (_gvcf_comp.is_range_compressable(known_pos_range2(si.pos,target_pos)))
        {
            assert(_block.count!=0);
            _block.count += (target_pos-_head_pos);
            _head_pos= target_pos;
        }
    }
}


void gvcf_aggregator::process(site_info& si)
{
    skip_to_pos(si.pos);
    add_site_internal(si);
}

void gvcf_aggregator::process(indel_info& ii)
{
    skip_to_pos(ii.pos);
    write_indel_record(ii);
    // TODO: eventually the block compression code would take this responsibility
    _last_indel.reset(new indel_info(ii));
}
void gvcf_aggregator::reset()
{
    _head.flush();
}

void gvcf_aggregator::flush()
{
    skip_to_pos(_report_range.end_pos);
    write_block_site_record();
}


//Add sites to queue for writing to gVCF
void
gvcf_aggregator::
add_site_internal(
    site_info& si)
{
    filter_site_by_last_indel_overlap(si);
    if (si.smod.is_phased_region)
    {
        _head_pos=si.pos+si.phased_ref.length();
    }
    else
    {
        _head_pos=si.pos+1;
    }
    // write_site
    queue_site_record(si);
}


void
gvcf_aggregator::
add_indel(const pos_t pos,
          const indel_key ik,
          const starling_diploid_indel_core& dindel,
          const starling_indel_report_info& iri,
          const starling_indel_sample_report_info& isri)
{
    indel_info ii(pos,ik,dindel,iri,isri);
    _head.process(ii);
}

// queue site record for writing, after
// possibly joining it into a compressed non-variant block
//
void
gvcf_aggregator::
queue_site_record(
    const site_info& si)
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
get_visible_alt_order(
    const site_info& si,
    std::vector<uint8_t>& altOrder)
{
    altOrder.clear();

    // list max_gt alts first:
    for (unsigned b(0); b<N_BASE; ++b)
    {
        if (b==si.dgt.ref_gt) continue;
        if (! DIGT::expect2(b,si.smod.max_gt)) continue;
        altOrder.push_back(b);
    }

#if 0
    // include other alts based on known count:
    for (unsigned b(0); b<N_BASE; ++b)
    {
        if (b==si.dgt.ref_gt) continue;
        if (DIGT::expect2(b,si.smod.max_gt)) continue;
        if (si.known_counts[b] > 0) altOrder.push_back(b);
    }
#endif
}



static
void
print_vcf_alt(
    const std::vector<uint8_t>& altOrder,
    std::ostream& os)
{
    bool is_print(false);
    for (const auto& b : altOrder)
    {
        if (is_print) os << ',';
        os << id_to_base(b);
        is_print=true;
    }
    if (! is_print) os << '.';
}



static
void
print_site_ad(
    const site_info& si,
    const std::vector<uint8_t>& altOrder,
    std::ostream& os)
{
    os << si.known_counts[si.dgt.ref_gt];

    for (const auto& b : altOrder)
    {
        os << ',' << si.known_counts[b];
    }
}


//writes out a SNP or block record
void
gvcf_aggregator::
write_site_record(
    const site_info& si) const
{
    std::ostream& os(*_osptr);

    os << _chrom << '\t'  // CHROM
       << (si.pos+1) << '\t'  // POS
       << ".\t";           // ID

    if (si.smod.is_phased_region)
    {
        os  << si.phased_ref << '\t'; // REF
    }
    else
    {
        os  << si.ref << '\t'; // REF
    }

    // ALT
    std::vector<uint8_t> altOrder;
    const bool isNoAlt(si.smod.is_unknown || (si.smod.is_block));
    if (isNoAlt)
    {
        os << '.';
    }
    else if (si.smod.is_phased_region)
    {
        os << si.phased_alt;
    }
    else
    {
        get_visible_alt_order(si,altOrder);
        print_vcf_alt(altOrder,os);
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
    if (! isNoAlt)
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
        if (si.smod.is_block)
        {
            os << _block.block_gqx.min();
        }
        else
        {
            if (si.smod.Qscore>=0)
            {
                os << si.smod.Qscore;
            }
            else
            {
                os << si.smod.gqx;
            }
        }
    }
    else
    {
        os << '.';
    }
    os << ':';
    //print DP:DPF
    if (si.smod.is_block)
    {
        os << _block.block_dpu.min() << ':'
           << _block.block_dpf.min();
    }
    else
    {
        os << si.n_used_calls << ':'
           << si.n_unused_calls;
    }

    if (isNoAlt) {}
    else if (si.smod.is_phased_region)
    {
        os << ':' << si.phased_AD;
    }
    else
    {
        os << ':';
        print_site_ad(si, altOrder, os);
    }
    os << '\n';
}






void
gvcf_aggregator::write_indel_record(const indel_info& ii)
{
    // flush any non-variant block before starting:
    write_block_site_record();

    std::ostream& os(*_osptr);

    os << _chrom << '\t'   // CHROM
       << ii.pos << '\t'   // POS
       << ".\t"            // ID
       << ii.iri().vcf_ref_seq << '\t'; // REF

    // ALT

    for (unsigned i = 0; i <ii._iri.size(); ++i)
    {
        if (i > 0) os << ',';
        os << ii.iri(i).vcf_indel_seq;
    }
    os << '\t';

    os << ii.dindel.indel_qphred << '\t'; //QUAL

    // FILTER:
    ii.imod().write_filters(os);
    os << '\t';

    // INFO
    os << "CIGAR=";
    for (unsigned i=0; i < ii._imod.size(); ++i)
    {
        if (i > 0) os << ',';
        os << ii.imod(i).cigar;
    }
    os << ';';
    os << "RU=";
    for (unsigned i = 0; i <ii._iri.size(); ++i)
    {
        if (i > 0) os << ',';
        if (ii.iri(i).is_repeat_unit() && ii.iri(i).repeat_unit.size() <= 20)
        {
            os << ii.iri(i).repeat_unit;
        }
        else
        {
            os << '.';
        }
    }
    os << ';';
    os << "REFREP=";
    for (unsigned i = 0; i <ii._iri.size(); ++i)
    {
        if (i > 0) os << ',';

        if (ii.iri(i).is_repeat_unit())
        {
            os << ii.iri(i).ref_repeat_count;
        }
        else
        {
            os << '.';
        }
    }
    os << ';';
    os << "IDREP=";
    for (unsigned i = 0; i <ii._iri.size(); ++i)
    {
        if (i > 0) os << ',';

        if (ii.iri(i).is_repeat_unit())
        {
            os << ii.iri(i).indel_repeat_count;
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
       << ii.imod().gq << ':';

    if (ii.imod().Qscore>=0)
        os << ii.imod().Qscore  << ':';
    else
        os << ii.imod().gqx  << ':';

    os << ii.isri().depth << ':';

    // SAMPLE AD:
    unsigned ref_count(0);
    for (unsigned i = 0; i <ii._isri.size(); ++i)
    {
        ref_count = std::max(ref_count, ii.isri(i).n_q30_ref_reads);
    }
    os << ref_count;
    for (unsigned i = 0; i <ii._isri.size(); ++i)
    {
        os << ',' << ii.isri(i).n_q30_indel_reads;
    }
    os << '\n';
}


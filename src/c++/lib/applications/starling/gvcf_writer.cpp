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

#include "gvcf_header.hh"
#include "blt_util/blt_exception.hh"
#include "blt_util/chrom_depth_map.hh"
#include "blt_util/io_util.hh"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include "gvcf_writer.hh"
#include "calibration_models.hh"
#include "indel_overlapper.hh"




//#define DEBUG_GVCF


#ifdef DEBUG_GVCF
#include "blt_util/log.hh"
#endif




void gvcf_writer::write_block_site_record()
{
    if (_block.count<=0) return;
    write_site_record(_block);
    _block.reset();
}

gvcf_writer::
gvcf_writer(
    const starling_options& opt,
    const starling_deriv_options& dopt,
    const reference_contig_segment& ref,
    const RegionTracker& nocompress_regions,
    std::ostream* osptr,
    const calibration_models& cm)
    : _opt(opt)
    , _report_range(dopt.report_range.begin_pos,dopt.report_range.end_pos)
    , _ref(ref)
    , _osptr(osptr)
    , _chrom(opt.bam_seq_name.c_str())
    , _dopt(dopt.gvcf)
    , _block(_opt.gvcf)
    , _head_pos(dopt.report_range.begin_pos)
    , _gvcf_comp(opt.gvcf,nocompress_regions)
    , _CM(cm)
{
    assert(_report_range.is_begin_pos);
    assert(_report_range.is_end_pos);

    if (! opt.gvcf.is_gvcf_output())
        throw std::invalid_argument("gvcf_writer cannot be constructed with nothing to do.");

    assert(nullptr != _osptr);
    assert((nullptr !=_chrom) && (strlen(_chrom)>0));

    if (! _opt.gvcf.is_skip_header)
    {
        finish_gvcf_header(_opt,_dopt, _dopt.chrom_depth,dopt.bam_header_data, *_osptr, cm);
    }

    variant_prefilter_stage::add_site_modifiers(_empty_site, _empty_site.smod, cm);

}


void gvcf_writer::filter_site_by_last_indel_overlap(digt_site_info& si)
{
    if (_last_indel)
    {
        if (si.pos >= _last_indel->end())
        {
            _last_indel.reset(nullptr);
        }
        else
        {
            indel_overlapper::modify_overlapping_site(*_last_indel, si, _CM);
        }
    }
}

// fill in missing sites
void
gvcf_writer::
skip_to_pos(const pos_t target_pos)
{
    // advance through any indel region by adding individual sites
    while (_head_pos<target_pos)
    {
        digt_site_info si = get_empty_site(_head_pos);

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


void gvcf_writer::process(std::unique_ptr<site_info> si)
{
    skip_to_pos(si->pos);

    if (typeid(*si) == typeid(digt_site_info))
    {
        add_site_internal(*downcast<digt_site_info>(std::move(si)));
    }
    else
    {
        add_site_internal(*downcast<continuous_site_info>(std::move(si)));
    }

}

void gvcf_writer::process(std::unique_ptr<indel_info> ii)
{
    skip_to_pos(ii->pos);

    if (typeid(*ii) == typeid(digt_indel_info))
    {
        auto ii_digt(downcast<digt_indel_info>(std::move(ii)));

        write_indel_record(*ii_digt);
        _last_indel = std::move(ii_digt);
    }
    else
    {
        write_indel_record(*downcast<continuous_indel_info>(std::move(ii)));
    }
}

void gvcf_writer::flush()
{
    skip_to_pos(_report_range.end_pos);
    write_block_site_record();
}


//Add sites to queue for writing to gVCF
void
gvcf_writer::
add_site_internal(
    digt_site_info& si)
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
gvcf_writer::
add_site_internal(
    continuous_site_info& si)
{
    // TODO: phasing
    _head_pos=si.pos+1;
    // write_site
    queue_site_record(si);
}




static
void
get_visible_alt_order(
    const digt_site_info& si,
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
    const digt_site_info& si,
    const std::vector<uint8_t>& altOrder,
    std::ostream& os)
{
    os << si.known_counts[si.dgt.ref_gt];

    for (const auto& b : altOrder)
    {
        os << ',' << si.known_counts[b];
    }
}

static
void
print_site_ad(
    const continuous_site_info& si,
    const continuous_site_call& call,
    std::ostream& os)
{
    os << si.known_counts[base_to_id(si.ref)] << "," << call._alleleDepth;
}



//writes out a SNP or block record
void
gvcf_writer::
write_site_record(
    const digt_site_info& si) const
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
    const bool isNoAlt(si.smod.is_unknown);
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
    if (si.dgt.is_snp)
    {
        os << "SNVSB=";
        {
            const StreamScoper ss(os);
            os << std::fixed << std::setprecision(1) << si.smod.strand_bias;
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
        if (si.smod.is_phasing_insufficient_depth)
        {
            os << ";Unphased";
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
    os << '\t';

    const bool is_nonref_gt(si.smod.max_gt != si.dgt.ref_gt);
    const bool is_print_pl(is_nonref_gt || si.dgt.is_snp || si.forcedOutput);

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
    if (is_print_pl)
    {
        os << ":PL";
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
        if (si.smod.Qscore>=0)
        {
            os << si.smod.Qscore;
        }
        else
        {
            os << si.smod.gqx;
        }
    }
    else
    {
        os << '.';
    }
    os << ':';
    //print DP:DPF
    os << si.n_used_calls << ':'
       << si.n_unused_calls;

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

    if (is_print_pl)
    {
        // print PL values
        os << ':';
        if (si.is_hetalt())
        {
            const unsigned print_gt(si.smod.max_gt);
            const uint8_t a0(DIGT::get_allele(print_gt,0));
            const uint8_t a1(DIGT::get_allele(print_gt,1));
            os << si.dgt.phredLoghood[si.dgt.ref_gt] << ','
               << si.dgt.phredLoghood[DIGT::get_gt_with_alleles(si.dgt.ref_gt,a0)] << ','
               << si.dgt.phredLoghood[DIGT::get_gt_with_alleles(a0,a0)] << ','
               << si.dgt.phredLoghood[DIGT::get_gt_with_alleles(si.dgt.ref_gt,a1)] << ','
               << si.dgt.phredLoghood[DIGT::get_gt_with_alleles(a0,a1)] << ','
               << si.dgt.phredLoghood[DIGT::get_gt_with_alleles(a1,a1)];
        }
        else if (si.dgt.is_haploid() || (si.smod.modified_gt == MODIFIED_SITE_GT::ONE))
        {
            os << si.dgt.phredLoghood[si.dgt.ref_gt] << ','
               << si.dgt.phredLoghood[si.smod.max_gt];
        }
        else
        {
            const unsigned print_gt(si.smod.max_gt);
            const uint8_t a0(DIGT::get_allele(print_gt,0));
            const uint8_t a1(DIGT::get_allele(print_gt,1));
            uint8_t alt(a0);
            if (si.dgt.ref_gt == a0)
            {
                alt = a1;
            }
            os << si.dgt.phredLoghood[si.dgt.ref_gt] << ','
               << si.dgt.phredLoghood[DIGT::get_gt_with_alleles(si.dgt.ref_gt,alt)] << ','
               << si.dgt.phredLoghood[DIGT::get_gt_with_alleles(alt,alt)];
        }
    }

    os << '\n';
}

void
gvcf_writer::write_site_record(
    const gvcf_block_site_record& si) const
{
    std::ostream& os(*_osptr);

    os << _chrom << '\t'  // CHROM
       << (si.pos+1) << '\t'  // POS
       << ".\t";           // ID

    os  << si.ref << '\t'; // REF

    // ALT
    os << '.';
    os << '\t';

    // QUAL:
    os << '.';
    os << '\t';

    // FILTER:
    si.write_filters(os);
    os << '\t';

    // INFO:
    if (si.count>1)
    {
        os << "END=" << (si.pos+si.count) << ';';
        os << _dopt.block_label;
    }
    else
    {
        os << '.';
    }
    os << '\t';

    //FORMAT
    os << "GT";
    os << ":GQX:DP:DPF";
    os << '\t';

    //SAMPLE
    os << si.get_gt() << ':';
    if (si.has_call)
    {
        os << _block.block_gqx.min();
    }
    else
    {
        os << '.';
    }
    os << ':';
    //print DP:DPF
    os << _block.block_dpu.min() << ':'
       << _block.block_dpf.min();
    os << '\n';
}



void
gvcf_writer::write_indel_record(const continuous_indel_info& ii)
{
    // flush any non-variant block before starting:
    write_block_site_record();

    std::ostream& os(*_osptr);


    for (auto& call : ii.calls)
    {

        os << _chrom << '\t'   // CHROM
           << ii.pos << '\t'   // POS
           << ".\t"            // ID
           << call._iri.vcf_ref_seq << '\t'; // REF

        // ALT
        os << call._iri.vcf_indel_seq;
        os << '\t';

        os << call.gq << '\t'; //QUAL

        // FILTER:
        call.write_filters(os);
        os << '\t';

        // INFO
        os << "CIGAR=";
        os << call.cigar;
        os << ';';
        os << "RU=";
        if (call._iri.is_repeat_unit() && call._iri.repeat_unit.size() <= 20)
        {
            os << call._iri.repeat_unit;
        }
        else
        {
            os << '.';
        }
        os << ';';
        os << "REFREP=";
        if (call._iri.is_repeat_unit())
        {
            os << call._iri.ref_repeat_count;
        }

        os << ';';
        os << "IDREP=";
        if (call._iri.is_repeat_unit())
        {
            os << call._iri.indel_repeat_count;
        }


        os << '\t';

        //FORMAT
        os << "GT:GQ:GQX:DPI:AD:VF" << '\t';

        //SAMPLE
        os << ii.get_gt() << ':'
           << call.gq << ':';

        os << call.gqx  << ':';

        os << call._isri.depth << ':';

        // SAMPLE AD:
        os << call._totalDepth - call._alleleDepth
           << ","
           << call._alleleDepth;
        // VF
        {
            const StreamScoper ss(os);
            os << ':' << std::setprecision(3) << call.variant_frequency();
        }
        os << '\n';
    }

}


void
gvcf_writer::write_indel_record(const digt_indel_info& ii)
{
    // flush any non-variant block before starting:
    write_block_site_record();

    std::ostream& os(*_osptr);
    auto& call(ii.first());

    os << _chrom << '\t'   // CHROM
       << ii.pos << '\t'   // POS
       << ".\t"            // ID
       << call._iri.vcf_ref_seq << '\t'; // REF

    // ALT

    for (unsigned i = 0; i <ii._calls.size(); ++i)
    {
        if (i > 0) os << ',';
        os << ii._calls[i]._iri.vcf_indel_seq;
    }
    os << '\t';

    os << call._dindel.indel_qphred << '\t'; //QUAL

    // FILTER:
    call.write_filters(os);
    os << '\t';

    // INFO
    os << "CIGAR=";
    for (unsigned i=0; i < ii._calls.size(); ++i)
    {
        if (i > 0) os << ',';
        os << ii._calls[i].cigar;
    }
    os << ';';
    os << "RU=";
    for (unsigned i = 0; i <ii._calls.size(); ++i)
    {
        const auto& iri(ii._calls[i]._iri);
        if (i > 0) os << ',';
        if (iri.is_repeat_unit() && iri.repeat_unit.size() <= 20)
        {
            os << iri.repeat_unit;
        }
        else
        {
            os << '.';
        }
    }
    os << ';';
    os << "REFREP=";
    for (unsigned i = 0; i <ii._calls.size(); ++i)
    {
        const auto& iri(ii._calls[i]._iri);
        if (i > 0) os << ',';

        if (iri.is_repeat_unit())
        {
            os << iri.ref_repeat_count;
        }
        else
        {
            os << '.';
        }
    }
    os << ';';
    os << "IDREP=";
    for (unsigned i = 0; i <ii._calls.size(); ++i)
    {
        const auto& iri(ii._calls[i]._iri);

        if (i > 0) os << ',';

        if (iri.is_repeat_unit())
        {
            os << iri.indel_repeat_count;
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
    os << "GT:GQ:GQX:DPI:AD:PL" << '\t';

    //SAMPLE
    os << ii.get_gt() << ':'
       << call.gq << ':';

    if (call.Qscore>=0)
        os << call.Qscore  << ':';
    else
        os << call.gqx  << ':';

    os << call._isri.depth << ':';

    // SAMPLE AD:
    unsigned ref_count(0);
    for (unsigned i = 0; i <ii._calls.size(); ++i)
    {
        ref_count = std::max(ref_count, ii._calls[i]._isri.n_q30_ref_reads);
    }
    os << ref_count;
    for (unsigned i = 0; i <ii._calls.size(); ++i)
    {
        os << ',' << ii._calls[i]._isri.n_q30_indel_reads;
    }

    // PL field
    os << ":";
    const unsigned icount(ii._calls.size());
    if (icount == 1)
    {
        using namespace STAR_DIINDEL;
        const auto& dindel(ii._calls[0]._dindel);
        const auto& pls(dindel.phredLoghood);
        if (dindel.is_haploid())
        {
            os << pls[NOINDEL] << ','
               << pls[HOM];
        }
        else
        {
            os << pls[NOINDEL] << ','
               << pls[HET] << ','
               << pls[HOM];
        }
    }
    else if (icount == 2)
    {
        // very roughly approximate the overlapping indel PL values
        //
        // 1. 0/0 - this is always maxQ
        // 2. 0/1 - set ot 0/0 from indel1
        // 3. 1/1 - set to 1/1 from indel0
        // 4. 0/2 - set to 0/0 from indel0
        // 5. 1/2 - this is always 0
        // 6. 2/2 - set to 1/1 from indel1
        //
        using namespace STAR_DIINDEL;
        const auto& pls0(ii._calls[0]._dindel.phredLoghood);
        const auto& pls1(ii._calls[1]._dindel.phredLoghood);

        os << starling_diploid_indel_core::maxQ << ','
           << pls1[NOINDEL] << ','
           << pls0[HOM] << ','
           << pls0[NOINDEL] << ','
           << 0 << ','
           << pls1[HOM];
    }
    else
    {
        assert(false && "Unexpected indel count");
    }

    os << '\n';
}

void
gvcf_writer::write_site_record(
    const continuous_site_info& si) const
{
    bool site_is_nonref = si.is_nonref();
    auto ref_base_id = base_to_id(si.ref);

    for (auto& call : si.calls)
    {
        bool is_no_alt(call._base == ref_base_id);

        // do not output the call for reference if the site has variants unless it is forced output
        if (!si.forcedOutput && site_is_nonref && is_no_alt)
            continue;

        std::ostream& os(*_osptr);

        os << _chrom << '\t'  // CHROM
           << (si.pos+1) << '\t'  // POS
           << ".\t";           // ID

        os  << si.ref << '\t'; // REF

        std::string gt(si.get_gt(call));

        // ALT
        if (is_no_alt)
            os << ".";
        else
            os << id_to_base(call._base);
        os << '\t';

        // QUAL: TODO - need separate calc?
        os << call.gq << '\t';

        // FILTER:
        call.write_filters(os);
        os << '\t';

        // INFO
        std::ostringstream info;

        if (si._is_snp)
        {
            info << "SNVSB=";
            {
                const StreamScoper ss(info);
                info << std::fixed << std::setprecision(1) << call.strand_bias;
            }
            info << ';';
            info << "SNVHPOL=" << si.hpol;
        }

        if (!is_no_alt)
        {
            if (_opt.do_codon_phasing)
            {
                if (!info.str().empty())
                    info << ";";
                info << "Unphased"; // TODO: placeholder until we do phasing on continuous variants
            }
        }

        os << (info.str().empty() ? "." : info.str()) << "\t";

        //FORMAT
        os << "GT";
        os << ":GQ";
        os << ":GQX";
        os << ":DP:DPF";
        if (!is_no_alt)
            os << ":AD";
        os << ":VF";

        os << '\t';

        //SAMPLE
        os << gt << ':';
        os << call.gq << ':';
        os << call.gqx << ':';
        // DP:DPF
        os << si.n_used_calls << ':' << si.n_unused_calls;

        if (!is_no_alt)
        {
            os << ':';
            print_site_ad(si, call, os);
        }

        {
            const StreamScoper ss(os);
            os << ':' << std::fixed << std::setprecision(3) << call.variant_frequency();
        }
        os << '\n';
    }
}

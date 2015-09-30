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

#include "gvcf_aggregator.hh"
#include "bedstreamprocessor.hh"
#include "variant_prefilter_stage.hh"
#include "indel_overlapper.hh"

gvcf_aggregator::
gvcf_aggregator(
    const starling_options& opt,
    const starling_deriv_options& dopt,
    const reference_contig_segment& ref,
    const RegionTracker& nocompress_regions,
    std::ostream* osptr,
    const pos_basecall_buffer& bc_buff)
    : _CM(opt, dopt.gvcf)
{
    if (! opt.gvcf.is_gvcf_output())
        throw std::invalid_argument("gvcf_aggregator cannot be constructed with nothing to do.");
    if (opt.is_ploidy_prior)
    {
        _writer.reset(new gvcf_writer(opt, dopt, ref, nocompress_regions, osptr, _CM));
        std::shared_ptr<variant_pipe_stage_base> overlapper(new indel_overlapper(_CM, ref, _writer));
        _codon_phaser.reset(new Codon_phaser(opt, bc_buff, ref, overlapper));
        std::shared_ptr<variant_pipe_stage_base> targeted_region_processor(new bed_stream_processor(opt.gvcf.targeted_regions_bedfile, opt.bam_seq_name.c_str(), _codon_phaser));
        _head.reset(new variant_prefilter_stage(_CM, targeted_region_processor));
    }
    else
    {
        // simpler pipeline when running in this mode
        _writer.reset(new gvcf_writer(opt, dopt, ref, nocompress_regions, osptr, _CM));
        std::shared_ptr<variant_pipe_stage_base> targeted_region_processor(new bed_stream_processor(opt.gvcf.targeted_regions_bedfile, opt.bam_seq_name.c_str(), _writer));
        _head.reset(new variant_prefilter_stage(_CM, targeted_region_processor));
    }
}

gvcf_aggregator::~gvcf_aggregator()
{
    _head->flush();
}

void
gvcf_aggregator::
add_site(std::unique_ptr<site_info> si)
{
    _head->process(std::move(si));
}

void
gvcf_aggregator::add_indel(std::unique_ptr<indel_info> info)
{
    _head->process(std::move(info));
}

void gvcf_aggregator::reset()
{
    _head->flush();
}





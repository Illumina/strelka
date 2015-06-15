
#include "gvcf_aggregator.hh"
gvcf_aggregator::
gvcf_aggregator(
    const starling_options& opt,
    const starling_deriv_options& dopt,
    const reference_contig_segment& ref,
    const RegionTracker& nocompress_regions,
    std::ostream* osptr,
    const pos_basecall_buffer& bc_buff)
    : _CM(opt, dopt.gvcf)
    , _writer(opt, dopt, ref, nocompress_regions, osptr, _CM)
    , _overlapper(_CM, ref, _writer)
    , _codon_phaser(opt, bc_buff, ref, _overlapper)
    , _targeted_region_processor(opt.gvcf.targeted_regions_bedfile, opt.bam_seq_name.c_str(), _codon_phaser)
    , _head(_CM, _targeted_region_processor)
{
    if (! opt.gvcf.is_gvcf_output())
        throw std::invalid_argument("gvcf_aggregator cannot be constructed with nothing to do.");


}

gvcf_aggregator::~gvcf_aggregator()
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

void gvcf_aggregator::reset()
{
    _head.flush();
}





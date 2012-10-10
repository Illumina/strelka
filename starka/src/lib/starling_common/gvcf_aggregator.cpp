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



gvcf_aggregator::
gvcf_aggregator(const starling_options& opt,
                const pos_range& report_range,
                std::ostream* osptr)
    : _opt(opt)
    , _report_range(report_range.begin_pos,report_range.end_pos)
    , _indel_end_pos(0)
    , _indel_buffer_size(0)
    , _osptr(osptr)
{
    assert(report_range.is_begin_pos);
    assert(report_range.is_end_pos);
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
            // buffer site:
            return;
        }
    }

    // write_site:
    *_osptr << _chrom << "\t" << (pos+1) << "\t" << "SNP" << "\n";
}



void
gvcf_aggregator::
add_indel(const pos_t pos,
          const starling_diploid_indel_core& dindel,
          const starling_indel_report_info& iri,
          const starling_indel_sample_report_info& isri) {


    // check to see if we add this indel to the buffer:
    if(0 != _indel_buffer_size) {
        // check if it's time to process the buffer:
        //if() {
            process_overlaps();
            //}
    }

    while(_indel_buffer.size() <= _indel_buffer_size) {
        _indel_buffer.push_back(indel_info());
    }
    _indel_buffer[_indel_buffer_size++].init(pos,dindel,iri,isri);
    // set indel end pos!!
}



void
gvcf_aggregator::
process_overlaps() {

    if(0==_indel_buffer_size) return;

    // do the overlap processing:

    _indel_buffer_size = 0;
}




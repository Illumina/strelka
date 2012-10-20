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

#ifndef __GVCF_AGGREGATOR_HH
#define __GVCF_AGGREGATOR_HH


#include "starling_common/gvcf_block_site_record.hh"
#include "starling_common/gvcf_locus_info.hh"

#include <iosfwd>

struct gvcf_deriv_options {
    
    gvcf_deriv_options()
        : is_max_depth(false)
        , max_depth(0)
    {}

    bool is_max_depth;
    double max_depth;
};



///
/// Assembles all site and indel call information into a consistent set, blocks output
/// and writes to a VCF stream
///

struct gvcf_aggregator {

    gvcf_aggregator(const starling_options& opt,
                    const pos_range& report_range,
                    const reference_contig_segment& ref,
                    std::ostream* os);

    ~gvcf_aggregator();

    void
    add_site(site_info& si);

    void
    add_indel(const pos_t pos,
              const indel_key ik,
              const starling_diploid_indel_core& dindel,
              const starling_indel_report_info& iri,
              const starling_indel_sample_report_info& isri);
    
    void
    flush() {
        skip_to_pos(_report_range.end_pos);
        process_overlaps();
        write_block_site_record();
    }

private:

    void
    add_site_internal(const site_info& si);

    void write_block_site_record() {
        if(_block.count<=0) return;
        write_site_record(_block.record);
        _block.reset();
    }

    void write_site_record(const site_info& si) const;

    void queue_site_record(const site_info& si);

    void modify_single_indel_record();

    void modify_overlap_indel_record();

    void modify_conflict_indel_record();

    // resolve a set of overlapping indel and site calls:
    void process_overlaps();
    
    void write_indel_record(const unsigned write_index=0);
    
    void
    skip_to_pos(const pos_t target_pos);

    const site_info&
    get_empty_site(const pos_t pos) {
        _empty_site.pos=pos;
        _empty_site.ref=_ref.get_base(pos);
        return _empty_site;
    }

    // initial policy is to write nothing at empty sites. why?
    //
    // (1) gatk does it
    // (2) if gVCF output is every turned off, the output would be ridiculous -- maybe this should be a gVCf only thing?
    //
    //void
    //write_empty_pos();

    //void
    //update_pos();

    const starling_options& _opt;
    const known_pos_range _report_range;
    const reference_contig_segment& _ref;
    std::ostream* _osptr;    // convenience to get gvcf stream from opt

    const char* _chrom;

    gvcf_deriv_options _dopt;

    pos_t _indel_end_pos;

    unsigned _indel_buffer_size;
    std::vector<indel_info> _indel_buffer;

    unsigned _site_buffer_size;
    std::vector<site_info> _site_buffer;

    gvcf_block_site_record _block;

    pos_t _head_pos; // we've observed sites upto but not includeing this position
    site_info _empty_site;
};


#endif

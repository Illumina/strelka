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


#include "blt_util/stream_stat.hh"
#include "starling_common/gvcf_locus_info.hh"

#include <iosfwd>



struct block_site_record {

    block_site_record(const starling_options& opt)
        : frac_tol(opt.gvcf_block_frac_tol)
        , abs_tol(opt.gvcf_block_abs_tol)
        , count(0)
    {}

    void
    reset() {
        count=0;
        block_gqx.reset();
    }

    // determine if the site could be joined this block:
    bool
    test(const site_info& si) const;

    // add record to this block
    void
    join(const site_info& si);

    site_info record;

    const double frac_tol;
    const int abs_tol;
    int count;
    stream_stat block_gqx;
    //stream_stat _blockDP;
    //stream_stat _blockMQ;
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
    set_chrom_name(const char* chrom) {
        _chrom=chrom;
    }

    void
    add_site(site_info& si);

    void
    add_indel(const pos_t pos,
              const indel_key ik,
              const starling_diploid_indel_core& dindel,
              const starling_indel_report_info& iri,
              const starling_indel_sample_report_info& isri);

private:

    void write_block_site_record() {
        if(_block.count<=0) return;
        write_site_record(_block.record);
        _block.reset();
    }

    void write_site_record(const site_info& si) const;

    void queue_site_record(site_info& si);

    void modify_single_indel_record();

    void modify_overlap_indel_record();

    void modify_conflict_indel_record();

    // resolve a set of overlapping indel and site calls:
    void process_overlaps();
    
    void write_indel_record(const unsigned write_index=0);


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

    pos_t _indel_end_pos;

    unsigned _indel_buffer_size;
    std::vector<indel_info> _indel_buffer;

    unsigned _site_buffer_size;
    std::vector<site_info> _site_buffer;

    block_site_record _block;
};


#endif

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


#include "starling_common/gvcf_locus_info.hh"

#include <iosfwd>


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
    add_site(const pos_t pos,
             const char ref,
             const unsigned n_used_calls,
             const unsigned n_unused_calls,
             const snp_pos_info& good_pi,
             const diploid_genotype& dgt,
             const bool is_nf_snp,
             const double sb,
             const unsigned hpol);

    void
    add_indel(const pos_t pos,
              const indel_key ik,
              const starling_diploid_indel_core& dindel,
              const starling_indel_report_info& iri,
              const starling_indel_sample_report_info& isri);

private:

    void modify_single_indel_record();

    void modify_overlap_indel_record();

    // resolve a set of overlapping indel and site calls:
    void process_overlaps();
    
    void write_indel_record();

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
    // convenient reference to gvcf stream from opt:
    std::ostream* _osptr;


    const char* _chrom;

    pos_t _indel_end_pos;

    unsigned _indel_buffer_size;
    std::vector<indel_info> _indel_buffer;

    unsigned _site_buffer_size;
    std::vector<site_info> _site_buffer;
};


#endif

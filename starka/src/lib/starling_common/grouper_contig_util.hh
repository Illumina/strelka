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

#ifndef __GROUPER_CONTIG_UTIL_HH
#define __GROUPER_CONTIG_UTIL_HH

#include "starling_common/grouper_contig.hh"

#include <iosfwd>
#include <memory>


/// returns false if genomic alignment is not possible:
///
bool
map_grouper_contig_read_to_genome(const grouper_contig& ctg,
                                  alignment& read_al);


bool
get_next_contig(std::istream& is,
                grouper_contig& ctg);


struct contig_reader {

    contig_reader(std::istream& is) : _is(is), _is_valid(true) {}

    bool next() {
        if(_is_valid) { _is_valid=get_next_contig(_is,_ctg); }
        return _is_valid;
    }
    
    const grouper_contig& get_contig() const { return _ctg; }

private:
    std::istream& _is;
    bool _is_valid;
    grouper_contig _ctg;
};




struct export_stream_reader;


// helps to manage io of contig and contig read file pair:
//
struct contig_data_manager {

    contig_data_manager(const std::string& contig_filename,
                        const std::string& contig_read_filename);

    ~contig_data_manager();

    export_stream_reader&
    contig_read_exr() { return *(_contig_read_exrp.get()); }

    contig_reader&
    creader() { return *(_creaderp.get()); }

private:
    std::auto_ptr<std::ifstream> _contig_isp;
    std::auto_ptr<std::ifstream> _contig_read_isp;
    std::auto_ptr<export_stream_reader> _contig_read_exrp;
    std::auto_ptr<contig_reader> _creaderp;
};


#endif

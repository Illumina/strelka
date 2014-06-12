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
/*
 * gvcfcompressor.h
 *
 *  Created on: Feb 21, 2014
 *      Author: Morten Kallberg
 */

#pragma once

#include <map>
#include "starling_common/gvcf_locus_info.hh"
#include "blt_common/blt_shared.hh"


typedef std::map<int, bool> posmap;
typedef std::map<std::string,posmap> chrposmap;
class gvcf_compressor {
public:
    gvcf_compressor();
    void read_bed(const std::string& input_file, const std::string& chrom);
    bool is_minor_allele_site(const int pos); //determine if the reference is the minor allele
    bool is_site_compressable(const gvcf_options& opt, const site_info& si); //determine if a site should be added to current block or if the block should be written out
    //determine the maximum no-call range that can be compressed, return the number of loci the block can be extended by
    int max_compressible_nocall_range(const int start, const int end);
    bool minor_allele_loaded;
//    void sanity_check();
private:
    chrposmap chr_to_pos;
    std::string my_chrom;
};

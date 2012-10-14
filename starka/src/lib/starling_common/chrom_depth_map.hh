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


#ifndef __CHROM_DEPTH_MAP
#define __CHROM_DEPTH_MAP


#include <map>
#include <string>


typedef std::map<std::string,double> cdmap_t;


// parse the chrom depth file
void
parse_chrom_depth(const std::string& chrom_depth_file,
                  cdmap_t& chrom_depth);


#endif

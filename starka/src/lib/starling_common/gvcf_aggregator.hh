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

///
/// Assembles all site and indel call information into a consistent set, blocks output
/// and writes to a VCF stream
///


// 1. naively cram everything together:
//
struct gvcf_aggregator {

    gvcf_aggregator(const starling_options& opt);

    add_site(site_info);

    add_indel(indel_info);

};


//std::ostream& operator<<(std::ostream& os,const alignment& al);


#endif

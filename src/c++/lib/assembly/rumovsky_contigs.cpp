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
 * assembler.cpp
 *
 *  Created on: Sep 10, 2014
 *  Author: Morten Kallberg
 */

#include <rumovsky_contigs.hh>
#include <array>
#include <sstream>
#include <vector>

#define DEBUG_RUM_CONTIG

#ifdef DEBUG_RUM_CONTIG
#include "blt_util/log.hh"
#endif

//rumovsky_contigs::rumovsky_contigs(){
//}
//
//Assembly rumovsky_contigs::do_assemble(){
//	Assembly as;
//	log_os << "Im assembling" << std::endl;
//
//	return as;
//}


// Add a SNP site to the assembly buffer


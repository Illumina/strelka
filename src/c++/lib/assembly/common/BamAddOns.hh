// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Rumovsky
// Copyright (c) 2013-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//
/// \file
/// \author Ole Schulz-Trieglaff
///

#pragma once

//#include <seqan/bam_io.h>
//#include <seqan/find.h>
//#include <seqan/seq_io.h>

#include "assembly/common/AssemblyReadInfo.hh"
#include "assembly/common/Qualities.hh"

// insertion operator for BamAlignment as Bamtools does not provide it
//std::ostream& operator<<(std::ostream& os, const seqan::BamAlignmentRecord& a);

//bool filterPoorAvgQuality(const seqan::BamAlignmentRecord& al, const int threshold);
//
//bool isAlignmentFiltered(const seqan::BamAlignmentRecord& align);
//
//bool isClippedAlignment(const seqan::BamAlignmentRecord& align);
//
//bool isIndelAlignment(const seqan::BamAlignmentRecord& align);
//
//bool isGoodShadow(const seqan::BamAlignmentRecord& align,
//                  const uint8_t lastMapq,
//                  const std::string& lastQname);

bool parseBamRegion(const std::string& bamRegion,
                    std::string& chrName,
                    unsigned& refStart,
                    unsigned& refEnd);

/**
 * Ensure that this interval on the reference has at least @p minLen
 */
void ensureRegionMinLength(const unsigned minLen, unsigned& start, unsigned& end);

/**
 * Loads the reference sequence from a file using region encoded as string
 */
bool loadRefSeq(const std::string& inRefFile,
                const std::string& bamRegion,
                std::string& refSeq /**< reference sequence is stored here */
               );

/**
 * Loads the reference sequence from a file using region given by chromosome, start and end
 */
bool loadRefSeq(const std::string& inRefFile,
                const std::string& chr,
                const int start,
                const int end,
                std::string& refSeq /**< reference sequence is stored here */
               );

//bool hasN(seqan::BamAlignmentRecord& bamAlign);
//
//bool isVariantRead(seqan::BamAlignmentRecord& record,
//                   uint8_t& lastMapq,
//                   std::string& lastQname);

bool readBamStream(const std::string& bamFile,
                   const std::string& bamRegion,
                   AssemblyReadInput& reads,
                   const bool onlyVariantReads);

void revComplement(std::string& fwd);

// quality trimming of read a la bwa
//void bwaQualityTrim(const int qualityCutoff, seqan::BamAlignmentRecord& record);

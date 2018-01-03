//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2018 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

#include "boost/test/unit_test.hpp"

#include "normalizeAlignment.hh"

#include "htsapi/align_path_bam_util.hh"
#include "htsapi/bam_util.hh"


//#define DEBUG_NORM_TEST


#ifdef DEBUG_NORM_TEST
#include "blt_util/log.hh"
#include <iostream>
#endif



BOOST_AUTO_TEST_SUITE( normalizeAlignment_test )


static
void
testNormalizeAlignment(
    const char* refSeq,
    const char* readSeq,
    const pos_t readPos,
    const char* cigar,
    const pos_t expectPos,
    const char* expectCigar)
{
#ifdef DEBUG_NORM_TEST
    log_os << "input ref: " << refSeq << " read: " << readSeq << " pos: " << readPos << " cigar: " << cigar << "\n";
#endif
    alignment al;
    al.pos=readPos;
    cigar_to_apath(cigar,al.path);
    string_bam_seq refBamSeq(refSeq,strlen(refSeq));
    string_bam_seq readBamSeq(readSeq,strlen(readSeq));

    const bool isChanged(normalizeAlignment(refBamSeq,readBamSeq,al));

#ifdef DEBUG_NORM_TEST
    log_os << "output pos: " << al.pos << " cigar: " << al.path << "\n";
#endif

    BOOST_REQUIRE(isChanged);
    BOOST_REQUIRE_EQUAL(al.pos, expectPos);
    BOOST_REQUIRE_EQUAL(apath_to_cigar(al.path),std::string(expectCigar));
}

BOOST_AUTO_TEST_CASE( test_normalizeAlignment )
{
    // basic left-shifting for hpols:
    testNormalizeAlignment("ACGTAAAATGTA","CGTAAATGT",1,"4M1D5M",1,"3M1D6M");
    testNormalizeAlignment("ACGTAAAATGTA","CGTAAAAATGT",1,"4M1I6M",1,"3M1I7M");
    testNormalizeAlignment("ACGTAAAATGTA","CGTAAAAAATGT",1,"3M1I2M1I5M",1,"3M2I7M");
    testNormalizeAlignment("ACGTAAAATGTA","CGTAATGT",1,"3M1D2M1D3M",1,"3M2D5M");

    // basic left-shifting for dinucs:
    testNormalizeAlignment("ACGTACACACACTGTA","CGTACACACTGT",1,"5M2D7M",1,"3M2D9M");
    testNormalizeAlignment("ACGTACACACACTGTA","CGTACACACACACTGT",1,"5M2I9M",1,"3M2I11M");

    // insert/delete collapsing
    testNormalizeAlignment("ACGTAAAATGTA","CGTAAAATGT",1,"4M1I1D5M",1,"10M");
    testNormalizeAlignment("ACGTAAAATGTA","CGTAAATGT",1,"4M1I2D4M",1,"3M1D6M");
    testNormalizeAlignment("ACGTAAAATGTA","CGTAAAAATGTA",1,"4M2I1D6M",1,"3M1I8M");

    // left-shift, then collapse:
    testNormalizeAlignment("ACGTAAAATGTA","CGTAAAAATG",1,"4M2I1M1D3M",1,"3M1I6M");
    // collapse, then left-shift
    testNormalizeAlignment("ATTGC","ATGC",0,"2M1I2D1M",0,"1M1D3M");

    // left-shift, merge indels, then more left-shifting for the merged indel
    //
    // This example is taken form HCC-1187BL hg19:chr2:27151600 and described in STREL-336,
    // in this case the second indel "8I" can be left-shifted adjacent to the first indel "6I".
    // When these indels are merged into one larger insertions "14I", that new merged insertion
    // can be left-shifted one more base. Failure to add the final left shift of the merged indel
    // was the cause of the error reported in STREL-336.
    testNormalizeAlignment("CACACACA","CACACACATATATACGCACACA",0,"3M6I3M8I2M",0,"2M14I6M");


    // The following two tests are motivated by a failure observed at hg19:chr2:179046342
    // in HCC1187BL leading to an unnormalized somatic indel, this is a second case addressed as
    // part of STREL-336. The identified issue was as follows: Whenever any indel was
    // left-shifted, the refPos and readPos values were not being updated to reflect the shift
    // and so the refPos/readPos pointers were effectively running ahead by the shift value.
    // If more than one indel is in the alignment, any left-shift computation on subsequent
    // indels was incorrect
    //

    // simplified demo based on the motivating issue
    testNormalizeAlignment("TAAAGT","TAAAAAAGAT",0,"4M3I1M1I1M",0,"1M3I4M1I1M");

    // motivating bug representing a failed read normalization aligned to hg19:chr2:179046342
    testNormalizeAlignment("AAAATATATATATATATATATA", "AATATATATATATATATATATATATATATATATATCTA", 0, "17M8I3M8I2M", 0,
                           "1M8I19M8I2M");


    // left-edge behavior
    testNormalizeAlignment("AAAAAAATGTA","AAAAAATGT",2,"3M1I5M",1,"9M");
    testNormalizeAlignment("AAAAAAATGTA","GGAAAAAATGT",2,"2S3M1I5M",1,"2S9M");
    testNormalizeAlignment("AAAAAATGTA","GGAAAATGT",1,"2S3M1D4M",2,"2S7M");
    testNormalizeAlignment("AAAAAATGTA","GGAAAATGT",2,"2S1D1I6M",2,"2S7M");
    testNormalizeAlignment("AAAAAAATGTA","AAAAAATGT",1,"1D1I8M",1,"9M");
    testNormalizeAlignment("AAAAAAATGTA","AAAAAATGT",1,"1H1S1D1I8M",1,"1H1S9M");

    // right-edge behavior
    testNormalizeAlignment("TGTGA","GTG",1,"3M1D",1,"3M");
    testNormalizeAlignment("TGTGA","GTGA",1,"3M1I",1,"4M");
    testNormalizeAlignment("TGTGA","GTGA",1,"3M1I1D",1,"4M");
    testNormalizeAlignment("TGTGA","GTGC",1,"3M1D1S",1,"3M1S");
    testNormalizeAlignment("TGTGA","GTGA",1,"3M1I1H",1,"4M1H");
    testNormalizeAlignment("TGTGA","GTGAC",1,"3M1I1D1S1H",1,"4M1S1H");

    // rna:
    testNormalizeAlignment("ACGTAAAATGTA","CGATGT",1,"2M3N1M1D4M",1,"2M3N1D5M");
    testNormalizeAlignment("ACGTAAAATGTA","CGAAATGT",1,"2M3N1M1I5M",1,"2M3N1I6M");

    // left-shifting policy that reduces mismatches:
    //
    // This example case suggested by Aaron Halpern in review.
    // The original left design would left shift only if the mismatch
    // count did not change, but it seems reasonable that left-shifting
    // if the mismatch count does not get worse is reasonable, and
    // consistent with the "collapsing of reference matching 1I1D patterns
    // as well (which is also a match-increasing transformation).
    //
    testNormalizeAlignment("TACAAAT","TAAAAT",0,"4M1D2M",0,"2M1D4M");
}



static
void
testNormalizeFromBamRecord(
    const char* refSeq,
    const char* readSeq,
    const pos_t readPos,
    const char* cigar,
    const pos_t expectPos,
    const char* expectCigar)
{
    reference_contig_segment ref;
    ref.seq() = refSeq;

    bam_record bamRecord;
    bam1_t& br(*bamRecord.get_data());
    bam1_core_t& ca(br.core);
    ca.pos = readPos;

    std::string bamRead(readSeq);
    uint8_t* qual = new uint8_t[bamRead.size()];
    edit_bam_read_and_quality(bamRead.c_str(),qual,br);
    delete[] qual;

    const std::string testCigar(cigar);
    ALIGNPATH::path_t inputPath;
    cigar_to_apath(testCigar.c_str(),inputPath);
    edit_bam_cigar(inputPath, br);

    const bool isChanged(normalizeBamRecordAlignment(ref,bamRecord));

    ALIGNPATH::path_t outputPath;
    bam_cigar_to_apath(bamRecord.raw_cigar(),bamRecord.n_cigar(),outputPath);

    BOOST_REQUIRE(isChanged);
    BOOST_REQUIRE_EQUAL(bamRecord.pos()-1,expectPos);
    BOOST_REQUIRE_EQUAL(apath_to_cigar(outputPath),std::string(expectCigar));
}

BOOST_AUTO_TEST_CASE( test_normalizeAlignmentFromBamRecord )
{
    testNormalizeFromBamRecord("ACGTAAAATGTA","CGTAAATGT",1,"4M1D5M",1,"3M1D6M");
    testNormalizeFromBamRecord("AAAAAAATGTA","AAAAAATGT",2,"3M1I5M",1,"9M");
}


BOOST_AUTO_TEST_SUITE_END()


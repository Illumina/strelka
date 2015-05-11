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


#include "BamAddOns.hh"

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

#include <seqan/seq_io.h>
#include <seqan/bam_io.h>
#include <seqan/modifier.h>

#include "common/Qualities.hh"

//static const uint16_t PHRED_OFFSET = 33;

// insertion operator for BamAlignment as Bamtools does not provide it
std::ostream& operator<<(std::ostream& os, const seqan::BamAlignmentRecord& a)
{
    os << a.qName << " " << a.seq << " " << a.rID << " " << a.beginPos << " " << a.mapQ;
    return os;
}

bool filterPoorAvgQuality(const seqan::BamAlignmentRecord& al, const int thr)
{
    // Hot fix for reads with empty quality strings, just wave them through
    if (length(al.qual) == 0) {
        return false;
    }
    const unsigned int readLen(length(al.seq));
    int avg(0);
    for (unsigned int i = 0; i < readLen; ++i)
    {
        avg += phred2char(al.qual[i]);
    }
    avg /= readLen;
    return (avg < thr);
}

bool isAlignmentFiltered(const seqan::BamAlignmentRecord& align)
{
    if ( hasFlagDuplicate(align) ||
         hasFlagQCNoPass(align)  ||
         hasFlagSecondary(align)
       )
    {
        return true;
    }
    return false;
}

bool isIndelAlignment(const seqan::BamAlignmentRecord& align)
{
    static const unsigned minIndelLen = 2;
    bool isBigIndel(false);

    for (auto cg : align.cigar)
    {
        if ( (cg.operation == 'D' && cg.count > minIndelLen) ||
             (cg.operation == 'I' && cg.count > minIndelLen)
           )
        {
            isBigIndel=true;
        }
    }
    return isBigIndel;
}

bool isClippedAlignment(const seqan::BamAlignmentRecord& align)
{

    static const unsigned minAvgQualClippedRead = 25;
    bool isGoodClipped(false);

    if (filterPoorAvgQuality(align,minAvgQualClippedRead))
    {
        return isGoodClipped;
    }

    const unsigned minClippedLen(5);
    for (auto cg : align.cigar)
    {
        //std::cout << cg.Type << " " << cg.Length << std::endl;
        if (cg.operation == 'S' && cg.count > minClippedLen)
        {
            isGoodClipped=true;
        }
    }
    return isGoodClipped;
}

bool isGoodShadow(const seqan::BamAlignmentRecord& align,
                  const uint8_t lastMapq,
                  const std::string& lastQname)
{

    // shadow read should be unmapped
    if (!hasFlagUnmapped(align) ) return false;
    // but its partner should be aligned
    if (hasFlagNextUnmapped(align)) return false;
    //if (!align.IsMateMapped()) return false;

    static const unsigned minAvgQualShadow = 25;
    static const uint16_t minSingletonMapq = 20;

    if (filterPoorAvgQuality(align,minAvgQualShadow))
    {
        return false;
    }

    if (align.qName != lastQname)
    {
        // something went wrong here, shadows should have their singleton partner
        // preceding them in the BAM file.
        return false;
    }

    if ((unsigned int)lastMapq > minSingletonMapq)
    {
        return true;
    }
    return false;
}

/**
 * Parses  the string representation of a bam region.
 */
bool parseBamRegion(const std::string& bamRegion,
                    std::string& chrName,
                    unsigned& start,
                    unsigned& end)
{

    boost::char_separator<char> sep(":-");
    boost::tokenizer< boost::char_separator<char> > tokTok(bamRegion, sep);
    std::vector<std::string> words(tokTok.begin(),tokTok.end());
    if (words.size()!=3 || words[0] == "")
    {
        std::cerr << "Unexpected format for bamRegion " << bamRegion << std::endl;
        return false;
    }
    chrName = words[0];
    try
    {
        start = boost::lexical_cast<unsigned>(words[1]);
        end   = boost::lexical_cast<unsigned>(words[2]);
    }
    catch (boost::bad_lexical_cast& e)
    {
        std::cerr << "Error parsing bam region " << bamRegion << std::endl;
        return false;
    }
    return true;
}

void ensureRegionMinLength(const unsigned wobble,
                           unsigned& refStart,
                           unsigned& refEnd)
{
    /*unsigned len = refEnd-refStart >= 0 ? (refEnd-refStart) : 0;;
    if (len < minLen)
    {
        int diff = minLen-len;
        refStart = (refStart-(diff/2) >= 0) ? refStart-(diff/2) : 0;
        refEnd += (diff/2);
    }*/
    
    refStart = (refStart>wobble) ? refStart-wobble : 0;
    refEnd += wobble;
    //std::cout << "Extending reference to " << refStart << " and " << refEnd << "\n";
}

bool loadRefSeq(const std::string& inRefFile,
                const std::string& chr,
                const int start,
                const int end,
                std::string& refSeq /**< reference sequence is stored here */
               )
{

    seqan::FaiIndex faiIndex;
    int res = read(faiIndex, inRefFile.c_str());
    if (res != 0)
    {
        std::cerr << "ERROR: Could not read FAI index " << inRefFile << ".fai\n";
        res = build(faiIndex, inRefFile.c_str());
        if (res != 0)
        {
            std::cerr << "ERROR: Could not build the FAI index, giving up.\n";
            return false;
        }
    }
    unsigned idx = 0;
    if (!getIdByName(faiIndex, chr, idx))
    {
        std::cerr << "ERROR: FAI index has no entry for chromosome " << chr << ".\n";
        return false;
    }
    seqan::CharString seqChr;
    if (readRegion(seqChr, faiIndex, idx, start, end) != 0)
    {
        std::cerr << "ERROR: Could not load region " << chr << ":" << start << "-" << end << ".\n";
        return false;
    }
    toUpper(seqChr);
    refSeq = toCString(seqChr);
    return true;
}


bool loadRefSeq(const std::string& inRefFile,
                const std::string& bamRegion,
                std::string& refSeq /**< reference sequence is stored here */
               )
{

    seqan::FaiIndex faiIndex;
    int res = read(faiIndex, inRefFile.c_str());
    if (res != 0)
    {
        std::cerr << "ERROR: Could not read FAI index " << inRefFile << ".fai\n";
        res = build(faiIndex, inRefFile.c_str());
        if (res != 0)
        {
            std::cerr << "ERROR: Could not build the FAI index, giving up.\n";
            return false;
        }
    }
    std::string chr;
    unsigned start, end;
    if ( !parseBamRegion(bamRegion,chr,start,end))
    {
        return false;
    }

    unsigned idx = 0;
    if (!getIdByName(faiIndex, chr, idx))
    {
        std::cerr << "ERROR: FAI index has no entry for chromosome " << chr << ".\n";
        return false;
    }
    seqan::CharString seqChr;
    if (readRegion(seqChr, faiIndex, idx, start, end) != 0)
    {
        std::cerr << "ERROR: Could not load region " << bamRegion << ".\n";
        return false;
    }
    toUpper(seqChr);
    refSeq = toCString(seqChr);
    return true;
}

/**
 * Tests the sequence of a BAM alignment for N characters.
 */
bool hasN(seqan::BamAlignmentRecord& bamAlign)
{
    seqan::CharString nChar("N");

    seqan::Finder<seqan::CharString> finder(bamAlign.seq);
    seqan::Pattern<seqan::CharString, seqan::Simple> pattern(nChar);
    if (find(finder, pattern))
    {
        return true;
    }
    return false;
}

/**
 * Tests if the BAM alignment is a variant breakpoint candidate.
 */
bool isVariantRead(seqan::BamAlignmentRecord& record,
                   uint8_t& lastMapq,
                   std::string& lastQname)
{

    bool isShadowKeeper(false);
    bool isClipKeeper(false);
    bool isIndelKeeper(false);

    //std::cout << "Testing " << record << "\n";
    if (isAlignmentFiltered(record)) return false;

    // skip reads contains N in their sequence
    if (hasN(record)) return false;

    if (isIndelAlignment(record))
    {
        isIndelKeeper = true;
    }

    if (isClippedAlignment(record))
    {
        isClipKeeper = true;
    }

    if (isGoodShadow(record,lastMapq,lastQname))
    {
        isShadowKeeper = true;
    }

    lastMapq = record.mapQ;
    assign(lastQname,record.qName);

    //std::cout << "result " << (isShadowKeeper || isClipKeeper || isIndelKeeper) << "\n";

    return (isShadowKeeper ||
            isClipKeeper   ||
            isIndelKeeper 
           );
}



/**
 * Opens a stream to read from BAM file.
 */
bool readBamStream(const std::string& bamFile,
                   const std::string& bamRegion,
                   AssemblyReadInput& reads,
                   const bool onlyVariantReads)
{

    seqan::Stream<seqan::Bgzf> bamStreamIn;
    if (!open(bamStreamIn, bamFile.c_str(),"r"))
    {
        std::cerr << "ERROR: Could not open " << bamFile << " for reading.\n";
        return false;
    }

    std::string chr;
    unsigned start, end;
    parseBamRegion(bamRegion,chr,start,end);
    //std::cout << "readBamStream chr=" << chr << " start=" << start << " end=" << end << "\n";

    // Read BAI index.
    seqan::BamIndex<seqan::Bai> baiIndex;
    // some optimistic assumptions about the name of the bam index
    std::string bamIndexFile = bamFile + ".bai";
    if (read(baiIndex, bamIndexFile.c_str()) != 0)
    {
        std::cerr << "ERROR: Could not read BAI index file " << bamIndexFile << "\n";
        return false;
    }

    // Setup name store, cache, and BAM I/O context.
    typedef seqan::StringSet<seqan::CharString> TNameStore;
    typedef seqan::NameStoreCache<TNameStore>   TNameStoreCache;
    typedef seqan::BamIOContext<TNameStore>     TBamIOContext;
    TNameStore      nameStore;
    TNameStoreCache nameStoreCache(nameStore);
    TBamIOContext   context(nameStore, nameStoreCache);

    // Read header.
    seqan::BamHeader header;
    if (readRecord(header, context, bamStreamIn, seqan::Bam()) != 0)
    {
        std::cerr << "ERROR: Could not read header from BAM file \n";
        return false;
    }

    // Translate from reference name to rID.
    int rID = 0;
    if (!getIdByName(nameStore, chr, rID, nameStoreCache))
    {
        std::cerr << "ERROR: Reference sequence named " << chr << " not known.\n";
        return false;
    }

    bool hasAlignments = false;
    if (!jumpToRegion(bamStreamIn, hasAlignments, context, rID, start, end, baiIndex))
    {
        std::cerr << "ERROR: Could not jump to " << bamRegion << "\n";
        return false;
    }
    if (!hasAlignments)
        return false;  // No alignments here.

    seqan::BamAlignmentRecord record;
    uint8_t lastMapq(0);
    std::string lastQname;

    while (!atEnd(bamStreamIn))
    {
        if (readRecord(record, context, bamStreamIn, seqan::Bam()) != 0)
        {
            std::cerr << "ERROR: Could not read record from BAM file.\n";
            return 1;
        }

        //std::cout << "Retrieved read " <<  record.qName << " " << record.seq << " "  << record.beginPos << "\n";
        // If we are on the next reference or at the end already then we stop.
        if (record.rID == -1 || record.rID > rID)
            break;
        
        // check if we're past the en
        // beginPos is int
        if (record.beginPos >= 0) {
            if ( (unsigned)record.beginPos >= end) 
                continue;
        } else {
            continue;
        }


        // If we are left of the selected position then we skip this record.
        if (record.beginPos+seqan::getAlignmentLengthInRef(record) < start)
            continue;

        if ( onlyVariantReads && ! isVariantRead(record, lastMapq, lastQname)) 
        {
            continue;
        } 

        std::string tmp;
        move(tmp, record.seq);
        reads.push_back(std::make_pair(record.beginPos,tmp));
    }
    return true;
}

void revComplement(std::string& seq)
{
    reverse(seq.begin(),seq.end());
    for (std::string::iterator p = seq.begin(); p != seq.end(); ++p) {
        switch (*p) {
        case 'A':
           *p = 'T';
           break;
        case 'C':
           *p = 'G';
           break;
        case 'G':
            *p = 'C';
            break;
        case 'T':
            *p = 'A';
            break;
        case 'N':
            *p = 'N';
            break;
        } // end of switch
    } // end of for
}

void bwaQualityTrim(const int qualityCutoff, seqan::BamAlignmentRecord& record) 
{
    int endpoint = 0; // not inclusive
    int max = 0;
    int i = length(record.seq) - 1;
    int terminalScore = char2phred(record.qual[i]);
    // Only perform soft-clipping if the last base has qual less than qualTrim
    if(terminalScore >= qualityCutoff)
        return;
    
    int subSum = 0;
    while(i >= 0)
    {
        int ps = char2phred(record.qual[i]);
        int score = qualityCutoff - ps;
        subSum += score;
        if(subSum > max)
        {
            max = subSum;
            endpoint = i;
        }
        --i;
    }
    // Clip the read
    record.seq  = prefix(record.seq, endpoint);
    record.qual = prefix(record.seq, endpoint);
}


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

///
/// \author Ole Schulz-Trieglaff
///


#include "GreedyAssembler.hh"

#include "graph/DeBruijnGraph.hh"
#include "graph/GraphWalker.hh"

#include "boost/foreach.hpp"
#include "boost/unordered_map.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/utility.hpp"

#include <cassert>
#include <vector>
#include <fstream>
#include <limits>

// compile with this macro to get verbose output:
//#define DEBUG_GREEDY_ASBL


static
void
computeRead2CtgMappings(const Assembly contigs,
                        const AssemblyReadInput& /*reads*/,
                        const unsigned readCount,
                        const unsigned wordLength,
                        const std::vector<str_uint_map_t>& readWordOffsets,
                        AssemblyReadOutput& readInfo
                       )
{
    unsigned cnt(0);
    for (Assembly::const_iterator ctgIter = contigs.begin(); ctgIter != contigs.end(); ++ctgIter)
    {
        const unsigned contigSize(ctgIter->seq.size());
        for (unsigned readIndex(0); readIndex<readCount; ++readIndex)
        {
            const str_uint_map_t& readWordOffset(readWordOffsets[readIndex]);
            AssemblyReadInfo& rinfo(readInfo[readIndex]);
            // store all reads sharing k-mers of the current word length with the contig
            for (unsigned j(0); j<=(contigSize-wordLength); ++j)
            {
                const std::string word(ctgIter->seq.substr(j,wordLength));
                str_uint_map_t::const_iterator findIter = readWordOffset.find(word);
                if (findIter != readWordOffset.end())
                {
                    /*std::cout << "Contig is " << ctgIter->seq << " cnt= " << cnt << "\n";
                    std::cout << "Read is " << reads[readIndex].second << " readIndex=" << readIndex << "\n";
                    std::cout << "Contig word " << word << " at " << j << " occurs in read " << readIndex << " at " << findIter->second << "\n";*/
                    rinfo.contigIds[cnt]     = true;
                    rinfo.contigOffsets[cnt] = (j - findIter->second);
                    //std::cout << "Setting read offset to " << rinfo.contigOffsets[cnt] << "\n";
                    break;
                }
            }
        }
        ++cnt;
    }

}

static
void
countUnalignedReads(const GreedyAssemblerOptions& opt,
                    const AssemblyReadInput& reads,
                    AssemblyReadOutput& readInfo,
                    Assembly& contigs,
                    Assembly& finalContigs,
                    //const std::string& refSeq,
                    std::vector<str_uint_map_t>& readWordOffsets,
                    unsigned& unusedReads,
                    const unsigned wordLength) {

    const unsigned readCount(reads.size());

    // finally -- set isUsed and decrement unusedReads
    for (Assembly::iterator ctgIter = contigs.begin(); ctgIter != contigs.end(); ++ctgIter)
    {
        //std::cout << "Checking ctg " << *ctgIter << "\n";
        //std::cout << "readCount = " << readCount << "\n";
        // we don't enforce that the sink k-mer is at the end of the contig sequence
        /*if (ctgIter->seq.find(sink) != std::string::npos)
        {
            ctgIter->hasSinkKmer = true;
        }*/

        const unsigned contigSize(ctgIter->seq.size());
        for (unsigned readIndex(0); readIndex<readCount; ++readIndex)
        {
            //std::cout << "readIndex = " << readIndex << "\n";
            /*const*/ str_uint_map_t& readWordOffset(readWordOffsets[readIndex]);
            AssemblyReadInfo& rinfo(readInfo[readIndex]);

            if (rinfo.isUsed) continue;

            // store all reads sharing k-mers of the current word length with the contig
            for (unsigned j(0); j<=(contigSize-wordLength); ++j)
            {
                const std::string word(ctgIter->seq.substr(j,wordLength));
                size_t pos = reads[readIndex].second.find(word);
                //str_uint_map_t::const_iterator findIter = readWordOffset.find(word);
                if (pos != std::string::npos)
                {
                    //std::cout << "Contig is " << ctgIter->seq << "\n";
                    //std::cout << "Read is " << reads[readIndex].second << "\n";
                    //std::cout << "Contig word " << word << " at " << j << " occurs in read " << readIndex << " at " << pos << "\n";
                    rinfo.isUsed                             = true;
                    rinfo.contigIds[finalContigs.size()]     = true;
                    rinfo.contigOffsets[finalContigs.size()] = pos;
                    //std::cout << "Setting read offset to " << rinfo.contigOffsets[finalContigs.size()] << "\n";

                    assert(unusedReads != 0);
                    --unusedReads;
                    ++ctgIter->seedReadCount;
                    readWordOffset[word] = j;
                    break;
                }
            }
        }
        ctgIter->avgCoverage /= ctgIter->seq.size();
        // throw away short stuff
        if (ctgIter->seq.length() <= opt.minContigLength) continue;
        finalContigs.push_back(*ctgIter);
    }
}

static
unsigned
selectMostFrequentKmer(const rumovsky::DeBruijnGraph& g, std::string& maxFreqKmer)
{
    unsigned int maxFreq(0);
    for(str_uint_map_t::const_iterator ct = g.getKmerFreqTable().begin(); ct != g.getKmerFreqTable().end(); ++ct) {

        if (ct->second > maxFreq) {
            maxFreqKmer = ct->first;
            maxFreq     = ct->second;
        }
    }
    return maxFreq;
}

static
bool
buildContigs(
    const GreedyAssemblerOptions& opt,
    const AssemblyReadInput& reads,
    AssemblyReadOutput& readInfo,
    const unsigned wordLength,
    Assembly& finalContigs,
    unsigned& unusedReads,
    const unsigned iteration
)
{
#ifdef DEBUG_GREEDY_ASBL
    std::cerr << __PRETTY_FUNCTION__ << " In GredyAssembler::buildContig. word length=" << wordLength << " number of reads " << reads.size() << "\n";
    for (unsigned readIndex(0); readIndex<reads.size(); ++readIndex)
    {
        std::cerr << reads[readIndex].first << "  " << reads[readIndex].second << " used=" << readInfo[readIndex].isUsed << "\n";
    }
#endif

    rumovsky::DeBruijnGraph g(wordLength);
    static const bool addReverse(true);
    if (!g.build(reads,readInfo,addReverse))
    {
        std::cerr << "Graph construction failed.\n";
        return false;
    }

    std::string mostFreqKmer;
    selectMostFrequentKmer(g,mostFreqKmer);
    std::cout << "Setting source to most frequent k-mer " << mostFreqKmer << "\n";

    std::string refSeq(""); // dummy reference
    std::string src(""); // dummy source node
    if (!opt.outFullGraphPrefix.empty())
    {
        std::string fullOutFileName = opt.outFullGraphPrefix + "_iter" + std::to_string(iteration);
        g.dumpToFile(fullOutFileName,refSeq,src);
    }

    if (opt.doPruning)
    {
        g.prune(opt.kmerPruningThreshold);
    }

    if (!opt.outGraphPrefix.empty()) {
        std::string outFileName = opt.outGraphPrefix + "_iter" + std::to_string(iteration);
        g.dumpToFile(outFileName,refSeq,src);
    }

//    std::cerr << "kmer freq table for DB graph class with size " <<  g.getKmerFreqTable().size() <<  "\n";
//    std::ofstream gHash("ghash.out");
//    for (str_uint_map_t::const_iterator ct = g.getKmerFreqTable().begin();ct!=g.getKmerFreqTable().end();++ct)
//    {
//        std::cout << ct->first << " " << ct->second << "\n";
//        //gHash << ct->first << " " << ct->second << "\n";
//    }
//    gHash.close();
//    std::cout << "Adjacency list for DB graph class: \n";
//    for (str_vec_map_t::const_iterator ct1 = g.getAdjList().begin(); ct1 != g.getAdjList().end();++ct1)
//    {
//        std::cout << ct1->first << " ";
//        for (AdjList::const_iterator ct2 = ct1->second.begin(); ct2 != ct1->second.end(); ++ct2 )
//        {
//            std::cout << ct2->first << " " << ct2->second << " ";
//        }
//        std::cout << "\n";
//    }

    std::cout << "Starting greedy traversal " << "\n";
    rumovsky::GraphWalkerGreedy<rumovsky::DeBruijnGraph> walker(g);

    Assembly contigs;
    walker.walk(mostFreqKmer,contigs);
    std::cout << "Graph traversal resulted in " << contigs.size() << " contigs.\n";

    const unsigned readCount(reads.size());
    // a set of read hashes; each read hash stores the starting positions of all kmers in the read
    std::vector<str_uint_map_t> readWordOffsets(readCount);

    countUnalignedReads(opt, reads, readInfo, contigs, finalContigs, readWordOffsets, unusedReads, wordLength);

    computeRead2CtgMappings(contigs,reads,readCount,wordLength,readWordOffsets,readInfo);
    // don't need this anymore:
    readWordOffsets.clear();

    if (finalContigs.size() == 0) return false;
    return true;
}



void
runGreedyAssembler(
    const GreedyAssemblerOptions& opt,
    const AssemblyReadInput& reads,
    AssemblyReadOutput& assembledReadInfo,
    Assembly& contigs)
{
#ifdef DEBUG_GREEDY_ASBL
    std::cerr << __PRETTY_FUNCTION__ << " Starting greedy assembly with " << reads.size() << " reads.\n";
#endif
    //assert(alphabet.size()>1);

    assembledReadInfo.clear();
    contigs.clear();

    assembledReadInfo.resize(reads.size());

    unsigned wordLength(opt.minWordLength);
    unsigned unusedReads(reads.size());

    for (unsigned iteration(0); iteration < opt.maxAssemblyIterations; ++iteration)
    {
        if (unusedReads < opt.minSeedReads) {
            std::cerr << "Not enough seed reads. Aborting. \n";
            return;
        }

        const unsigned lastUnusedReads(unusedReads);
        for (; wordLength<=opt.maxWordLength; wordLength+=opt.wordStepSize)
        {
            if (opt.verbose)
            {
                std::cerr << "Starting assembly with k=" << wordLength << " iter=" << iteration << "\n";
                std::cerr << "opt.minWordLength=" << opt.minWordLength << " opt.maxWordLength " << opt.maxWordLength << "\n";
            }
            const bool isAssemblySuccess = buildContigs(opt, reads, assembledReadInfo, wordLength,
                                                        contigs, unusedReads, iteration);
            if (opt.verbose)
            {
                std::cerr << "isAssemblySuccess=" << isAssemblySuccess << "\n";
            }
            if (isAssemblySuccess) break;
        }

#ifdef DEBUG_GREEDY_ASBL
        std::cerr << __PRETTY_FUNCTION__ << " iter: " << iteration << " unused readMap now: " << unusedReads << "\n";
#endif

        // stop if no change in number of unused reads
        if (unusedReads == lastUnusedReads)
        {
#ifdef DEBUG_GREEDY_ASBL
            std::cerr << __PRETTY_FUNCTION__ << " Number of unused reads (" << unusedReads << ") did not change in this iteration. Stopping.\n";
#endif
            break;
        }

    }
/*#ifdef DEBUG_GREEDY_ASBL
    std::cerr << __PRETTY_FUNCTION__ << " Reached max number of assembly iterations: " << opt.maxAssemblyIterations << "\n";
#endif*/

}


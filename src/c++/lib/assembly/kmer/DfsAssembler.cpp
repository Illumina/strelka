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


#include "DfsAssembler.hh"

#include "assembly/graph/DeBruijnGraph.hh"
#include "assembly/graph/GraphWalker.hh"

#include "boost/foreach.hpp"
#include "boost/unordered_map.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/utility.hpp"

#include <cassert>
#include <vector>
#include <fstream>
#include <limits>

// compile with this macro to get verbose output:
//#define DEBUG_ASBL


//static
//bool
//setSourceInGraph(const std::string& refSeq,
//                 const rumovsky::DeBruijnGraph& g,
//                 std::string& source)
//{
//    unsigned const refSeqLen(refSeq.size());
//    const unsigned int wordLength(g.getWordLength());
//    assert(wordLength < refSeq.size());
//
//    //unsigned sourceRefIdx(0);
//    // set source
//    for (unsigned i(0); i<=(refSeqLen-wordLength); ++i)
//    {
//        std::string refKmer(refSeq.substr(i,wordLength));
//        if (g.isInGraph(refKmer))
//        {
//            source = refKmer;
//            //sourceRefIdx = i;
//            break;
//        }
//    }
//    if (source.empty())
//    {
//        std::cerr << " Empty source k-mer at k=" << wordLength  << ". Have " << g.size()  << " words left. \n";
//        return false;
//    }
//    std::cerr << "Setting graph source k-mer to : " << source << "\n";
//    return true;
//}
//
//
//// check if reference and contig have a kmer in common
//// quadratic runtime, could be improved
//static
//bool
//hasRefKmer(const std::string& refSeq,
//           const std::string& contig,
//           const unsigned wordLength)
//{
//    const unsigned refLen(refSeq.size());
//    const unsigned ctgLen(contig.size());
//
//    assert(ctgLen>=wordLength && refLen>=wordLength);
//
//    bool refKmerFound(false);
//
//    for (unsigned i=0;i<(refLen-wordLength);++i)
//    {
//        const std::string kmer = refSeq.substr(i,wordLength);
//        if (contig.find(kmer) != std::string::npos)
//        {
//            refKmerFound=true;
//            break;
//        }
//    }
//    return refKmerFound;
//}
//
//// check if kmer is contained in reference
///*static
//bool
//isRefKmer(const std::string& refSeq,
//           const std::string& kmer)
//{
//    return ( refSeq.find(kmer) != std::string::npos);
//}
//*/
//
//
//static
//void
//computeRead2CtgMappings(const Assembly contigs,
//                        const AssemblyReadInput& /*reads*/,
//                        const unsigned readCount,
//                        const unsigned wordLength,
//                        const std::vector<str_uint_map_t>& readWordOffsets,
//                        AssemblyReadOutput& readInfo
//                       )
//{
//    unsigned cnt(0);
//    for (Assembly::const_iterator ctgIter = contigs.begin(); ctgIter != contigs.end(); ++ctgIter)
//    {
//        const unsigned contigSize(ctgIter->seq.size());
//        for (unsigned readIndex(0); readIndex<readCount; ++readIndex)
//        {
//            const str_uint_map_t& readWordOffset(readWordOffsets[readIndex]);
//            AssemblyReadInfo& rinfo(readInfo[readIndex]);
//            // store all reads sharing k-mers of the current word length with the contig
//            for (unsigned j(0); j<=(contigSize-wordLength); ++j)
//            {
//                const std::string word(ctgIter->seq.substr(j,wordLength));
//                str_uint_map_t::const_iterator findIter = readWordOffset.find(word);
//                if (findIter != readWordOffset.end())
//                {
//                    /*std::cout << "Contig is " << ctgIter->seq << " cnt= " << cnt << "\n";
//                    std::cout << "Read is " << reads[readIndex].second << " readIndex=" << readIndex << "\n";
//                    std::cout << "Contig word " << word << " at " << j << " occurs in read " << readIndex << " at " << findIter->second << "\n";*/
//                    rinfo.contigIds[cnt]     = true;
//                    rinfo.contigOffsets[cnt] = (j - findIter->second);
//                    //std::cout << "Setting read offset to " << rinfo.contigOffsets[cnt] << "\n";
//                    break;
//                }
//            }
//        }
//        ++cnt;
//    }
//
//}
//
//static
//void
//countUnalignedReads(const DfsAssemblerOptions& opt,
//                    const AssemblyReadInput& reads,
//                    AssemblyReadOutput& readInfo,
//                    Assembly& contigs,
//                    Assembly& finalContigs,
//                    const std::string& refSeq,
//                    std::vector<str_uint_map_t>& readWordOffsets,
//                    unsigned& unusedReads,
//                    const unsigned wordLength) {
//
//    const unsigned readCount(reads.size());
//
//    // finally -- set isUsed and decrement unusedReads
//    for (Assembly::iterator ctgIter = contigs.begin(); ctgIter != contigs.end(); ++ctgIter)
//    {
//        //std::cout << "Checking ctg " << *ctgIter << "\n";
//        //std::cout << "readCount = " << readCount << "\n";
//        // we don't enforce that the sink k-mer is at the end of the contig sequence
//        if (hasRefKmer(refSeq,ctgIter->seq,wordLength))
//        {
//            ctgIter->hasSinkKmer = true;
//        }
//        /*if (ctgIter->seq.find(sink) != std::string::npos)
//        {
//            ctgIter->hasSinkKmer = true;
//        }*/
//
//        const unsigned contigSize(ctgIter->seq.size());
//        for (unsigned readIndex(0); readIndex<readCount; ++readIndex)
//        {
//            //std::cout << "readIndex = " << readIndex << "\n";
//            /*const*/ str_uint_map_t& readWordOffset(readWordOffsets[readIndex]);
//            AssemblyReadInfo& rinfo(readInfo[readIndex]);
//
//            if (rinfo.isUsed) continue;
//
//            // store all reads sharing k-mers of the current word length with the contig
//            for (unsigned j(0); j<=(contigSize-wordLength); ++j)
//            {
//                const std::string word(ctgIter->seq.substr(j,wordLength));
//                size_t pos = reads[readIndex].second.find(word);
//                //str_uint_map_t::const_iterator findIter = readWordOffset.find(word);
//                if (pos != std::string::npos)
//                {
//                    //std::cout << "Contig is " << ctgIter->seq << "\n";
//                    //std::cout << "Read is " << reads[readIndex].second << "\n";
//                    //std::cout << "Contig word " << word << " at " << j << " occurs in read " << readIndex << " at " << pos << "\n";
//                    rinfo.isUsed                             = true;
//                    rinfo.contigIds[finalContigs.size()]     = true;
//                    rinfo.contigOffsets[finalContigs.size()] = pos;
//                    //std::cout << "Setting read offset to " << rinfo.contigOffsets[finalContigs.size()] << "\n";
//
//                    assert(unusedReads != 0);
//                    --unusedReads;
//                    ++ctgIter->seedReadCount;
//                    readWordOffset[word] = j;
//                    break;
//                }
//            }
//        }
//        ctgIter->avgCoverage /= ctgIter->seq.size();
//        // throw away short stuff
//        if (ctgIter->seq.length() <= opt.minContigLength) continue;
//        finalContigs.push_back(*ctgIter);
//    }
//}
//
//static
//bool
//buildContigs(
//    const DfsAssemblerOptions& opt,
//    const AssemblyReadInput& reads,
//    AssemblyReadOutput& readInfo,
//    const unsigned wordLength,
//    Assembly& finalContigs,
//    unsigned& unusedReads,
//    const std::string& refSeq,
//    const unsigned iteration
//)
//{
//#ifdef DEBUG_ASBL
//    std::cerr << __PRETTY_FUNCTION__ << " In SmallAssembler::buildContig. word length=" << wordLength << " number of reads " << reads.size() << "\n";
//    for (unsigned readIndex(0); readIndex<reads.size(); ++readIndex)
//    {
//        std::cerr << reads[readIndex].first << "  " << reads[readIndex].second << " used=" << readInfo[readIndex].isUsed << "\n";
//    }
//#endif
//
//    rumovsky::DeBruijnGraph g(wordLength);
//    if (!g.build(reads,readInfo))
//    {
//        std::cerr << "Graph construction failed.\n";
//        return false;
//    }
//
//    // prune low-frequency kmers
//    if (opt.doPruning)
//    {
//        std::cerr << "Pruning graph with threshold " << opt.kmerPruningThreshold << "\n";
//        g.prune(opt.kmerPruningThreshold);
//    }
//
//    // set source node
//    std::string srcGraph;
//    if (!setSourceInGraph(refSeq,g,srcGraph))
//    {
//        return false;
//    }
//
//    // write pruned graph but before path contraction
//    if (!opt.outFullGraphPrefix.empty())
//    {
//        std::string outFullFileName = opt.outFullGraphPrefix + "_iter" + std::to_string(iteration);
//        g.dumpToFile(outFullFileName,refSeq,srcGraph);
//    }
//
//    // contract paths
//    g.contractPaths(srcGraph);
//    if (!opt.outGraphPrefix.empty())
//    {
//        std::string outFileName = opt.outGraphPrefix + "_iter" + std::to_string(iteration);
//        g.dumpToFile(outFileName,refSeq,srcGraph);
//    }
//    rumovsky::GraphWalkerDFS<rumovsky::DeBruijnGraph> walker(g);
//
//    Assembly contigs;
//    walker.walk(srcGraph,contigs);
//    std::cout << "Graph traversal resulted in " << contigs.size() << " contigs.\n";
//
//    /*std::cerr << "kmer freq table for DB graph class \n";
//    std::ofstream gHash("ghash.out");
//    for (str_uint_map_t::const_iterator ct = g.getKmerFreqTable().begin();ct!=g.getKmerFreqTable().end();++ct)
//    {
//        gHash << ct->first << " " << ct->second << "\n";
//    }
//    gHash.close();*/
//    /*std::cout << "Adjacency list for DB graph class: \n";
//    for (str_vec_map_t::const_iterator ct1 = g.getAdjList().begin(); ct1 != g.getAdjList().end();++ct1)
//    {
//        std::cout << ct1->first << " " << ct1->second.size() << " ";
//        for (AdjList::const_iterator ct2 = ct1->second.begin(); ct2 != ct1->second.end(); ++ct2 )
//        {
//            std::cout << ct2->first << " ";
//        }
//        std::cout << "\n";
//    }*/
//
//    const unsigned readCount(reads.size());
//    // a set of read hashes; each read hash stores the starting positions of all kmers in the read
//    std::vector<str_uint_map_t> readWordOffsets(readCount);
//
//    countUnalignedReads(opt, reads, readInfo, contigs, finalContigs, refSeq, readWordOffsets, unusedReads, wordLength);
//
//    computeRead2CtgMappings(contigs,reads,readCount,wordLength,readWordOffsets,readInfo);
//    // don't need this anymore:
//    readWordOffsets.clear();
//
//    if (finalContigs.size() == 0) return false;
//    return true;
//}

void
runDfsAssembler(
    const DfsAssemblerOptions& opt,
    const AssemblyReadInput& reads,
    const std::string& refSeq,
    AssemblyReadOutput& assembledReadInfo,
    Assembly& contigs)

{
//    assert(alphabet.size()>1);

    assembledReadInfo.clear();
    contigs.clear();

    assembledReadInfo.resize(reads.size());

    unsigned wordLength(opt.minWordLength);
    unsigned unusedReads(reads.size());

    for (unsigned iteration(0); iteration < opt.maxAssemblyIterations; ++iteration)
    {
        if (unusedReads < opt.minSeedReads) return;

        const unsigned lastUnusedReads(unusedReads);
        for (; wordLength<=opt.maxWordLength; wordLength+=opt.wordStepSize)
        {
            if (opt.verbose)
            {
                std::cerr << "Starting assembly with k=" << wordLength << " iter=" << iteration << "\n";
                std::cerr << "opt.minWordLength=" << opt.minWordLength << " opt.maxWordLength " << opt.maxWordLength << std::endl;
            }
            const bool isAssemblySuccess = true;
//            const bool isAssemblySuccess = true; buildContigs(opt, reads, assembledReadInfo, wordLength,
//                                                        contigs, unusedReads, refSeq, iteration);
            if (opt.verbose)
            {
                std::cerr << "isAssemblySuccess=" << isAssemblySuccess << std::endl;
            }
            if (isAssemblySuccess) break;
        }

        // stop if no change in number of unused reads
        if (unusedReads == lastUnusedReads)
        {
            break;
        }

    }

}


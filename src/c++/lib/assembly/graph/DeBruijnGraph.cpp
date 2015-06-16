#include "assembly/graph/DeBruijnGraph.hh"

#include <algorithm>  
#include <iterator>  

namespace rumovsky {

static
void
generateAllKmers(const std::string& str,
                 const unsigned k,
                 std::vector<std::string>& kmers
                 )
{
    static const unsigned len(str.size());
    assert(k<len);
    for (unsigned i(0);i<(len-k);++i) {
        kmers.push_back(str.substr(i,k));
    }
}

static
std::string
getAfix(const std::string& contig,
        const unsigned length,
        const AFFIX_TYPE af)
{

    const unsigned csize(contig.size());
    /*if(length <= csize) {
        std::cerr << "contig=" << contig << std::endl;
    }*/
    assert(length <= csize);

    if (af==AFFIX_TYPE::SUFFIX) return contig.substr((csize-length),length);
    else       return contig.substr(0,length);
}


static
std::string
addBase(const std::string& contig,
        const char base,
        const bool atEnd)
{
    return atEnd ? (contig + base) : (base + contig);
}

template <typename Iter, typename Cont>
bool is_last(Iter iter, const Cont& cont)
{
    return (iter != cont.end()) && (boost::next(iter) == cont.end());
}


bool DeBruijnGraph::build(const AssemblyReadInput& reads,
                          const AssemblyReadOutput& readInfo,
						  const bool addreverse)
{
    const unsigned readCount(reads.size());
    readWordOffsets_.resize(readCount);

//    for (unsigned readIndex(0); readIndex<readCount; ++readIndex)
//    {
//        const AssemblyReadInfo& rinfo(readInfo[readIndex]);
//        // skip reads used in a previous iteration
//        if (rinfo.isUsed) continue;
//
//        // stores the index of a kmer in a read sequence
//        const std::string& seq(reads[readIndex].second);
//        const unsigned readLen(seq.size());
//
//        // just in case
//        if (readLen < wordLength_) continue;
//
//        str_uint_map_t& readWordOffset(readWordOffsets_[readIndex]);
//
//        for (unsigned j(0); j<=(readLen-wordLength_); ++j)
//        {
//            std::string kmer(seq.substr(j,wordLength_));
//            if (readWordOffset.find(kmer) != readWordOffset.end())
//            {
//                // try again with different k-mer size
////#ifdef DEBUG_DB_GRAPH
//                std::cerr << "word " << kmer << " repeated in read " << readIndex << "\n";
//                //std::cerr << "res = " << readWordOffset.find(word)->second << "\n";
////#endif
//                return false;
//            }
//            // record (0-indexed) start point for word in read
//            #ifdef DEBUG_DB_GRAPH
//            std::cerr << "Recording " << kmer << " at " << j << " in " << seq <<  "\n";
//            #endif
//            //std::cerr << "Recording " << word << " at " << j << " in " << seq <<  "\n";
//            readWordOffset[kmer]=j;
//            // count occurrences
//            ++kmerFreqTable_[kmer];
//
//            if (addReverse) {
//                revComplement(kmer); // seqan
//                readWordOffset[kmer]=j;
//                // count occurrences
//                ++kmerFreqTable_[kmer];
//
//            }
//        }
//    }
//
//    // do this at the end
//    buildAdjList_();
    return true;
}

void DeBruijnGraph::prune(const unsigned minFreq)
{
    // prune hash, remove singleton k-mers
    unsigned pruneCnt(0);
    for (str_uint_map_t::iterator it=kmerFreqTable_.begin();
         it!=kmerFreqTable_.end();)
    {
        if (it->second <= minFreq)
        {
            const std::string wordToDelete(it->first);
            it = kmerFreqTable_.erase(it);
            ++pruneCnt;
        }
        else
        {
            ++it;
        }
    }
    // rebuild adjList from pruning k-mer table
    buildAdjList_();

//#ifdef DEBUG_DB_GRAPH
    std::cerr << "Pruned " << pruneCnt << " k-mers. Remaining " << kmerFreqTable_.size() << std::endl;
//#endif
}

void DeBruijnGraph::getNeighbours(const std::string& vertex, AdjList& adj) const
{
    str_vec_map_t::const_iterator cter = adjList_.find(vertex);
    if (cter != adjList_.end()){
        adj = cter->second;
    }
}

void DeBruijnGraph::buildAdjList_()
{
    adjList_.clear();

    for (str_uint_map_t::const_iterator ct = kmerFreqTable_.begin();
         ct != kmerFreqTable_.end(); ++ct)
    {

        for ( int afInt = PREFIX; afInt <= SUFFIX; ++afInt )
        {
            AFFIX_TYPE af = static_cast<AFFIX_TYPE>(afInt);
            std::string tmp(getAfix(ct->first,wordLength_-1,af));

            bool atEnd(af==AFFIX_TYPE::SUFFIX);
            static const std::string alphabet("ACGT");
            for(const auto c : alphabet)
            {
                const std::string overlap(addBase(tmp,c,atEnd));
                if (overlap == ct->first) continue; // prevent self-loops
                if (kmerFreqTable_.count(overlap))
                {
                    adjList_[ct->first].push_back(std::make_pair(overlap,af));
                    //std::cerr << "Adjacency : " << ct->first << " => " << overlap << "\n";
                }
            }
        }
        //std::cerr << "Node " << ct->first << " has " << adjList_[ct->first].size() << " edges. \n";
    }

}


void
DeBruijnGraph::contractPaths(const std::string& startVertex)
{
#ifdef DEBUG_DB_GRAPH
    std::cerr << "INIT MAX PATHS SEARCH START " << startVertex << "\n";
#endif
    str_bool_map_t seenVertices;
    seenVertices[startVertex] = true;
    Contig ctg;
    ctg.seq         = startVertex;
    ctg.avgCoverage = 0;
    ctg.numKmers    = 1;

    str_uint_map_t newKmerFreqTable;
    str_vec_map_t newAdjList;

    doPathContract_(ctg, seenVertices,newKmerFreqTable,newAdjList);

#ifdef DEBUG_DB_GRAPH
    std::cerr << "seenVertices=" << seenVertices.size() << " newKmerFreqTable=" <<  newKmerFreqTable.size() << "\n";
#endif

    kmerFreqTable_ = newKmerFreqTable;
    adjList_       = newAdjList;
}

std::string
DeBruijnGraph::doPathContract_(Contig& contigSoFar,
                               str_bool_map_t& seenVertices,
                               str_uint_map_t& newKmerFreqTable,
                               str_vec_map_t& newAdjList)
{

#ifdef DEBUG_DB_GRAPH
    std::cerr << "findMaxPaths : contigSoFar = " << contigSoFar.seq << "\n";
#endif
    assert(contigSoFar.seq.size()>=wordLength_);

    std::string lastKmer(getAfix(contigSoFar.seq,wordLength_,AFFIX_TYPE::SUFFIX));
    str_vec_map_t::const_iterator findIter = adjList_.find(lastKmer);

    if (findIter == adjList_.end())
    {
        std::cerr << "kmer " << lastKmer << " not in hash. \n";
        assert(0);
    }

    // Ideally this would happen in the API
    AdjList nbsFiltered;
    std::copy_if(findIter->second.begin(),findIter->second.end(), std::back_inserter(nbsFiltered),
                [](Neighbour const & n) { return n.second == AFFIX_TYPE::SUFFIX; } );

    unsigned numNeighbours = nbsFiltered.size();
    //assert(numNeighbours != 0);

    if (numNeighbours == 0) 
    {

#ifdef DEBUG_DB_GRAPH
        std::cerr << "Dead end, adding new contig " << contigSoFar.seq << std::endl;
#endif
        newKmerFreqTable[contigSoFar.seq] = contigSoFar.numKmers;
    }
    else if (numNeighbours == 1)
    {
        // in- and out-degree = 1
        bool movedForward(false);
        for (auto s: findIter->second)
        {
            // don't walk back
            if (seenVertices[s.first]) continue;
            if (s.second != AFFIX_TYPE::SUFFIX) continue;
            std::string lastNucleotide(getAfix(s.first,1,AFFIX_TYPE::SUFFIX));
#ifdef DEBUG_DB_GRAPH
            std::cerr << " lastKmer = " << lastKmer << " found neighbour = " << s.first << std::endl;
            std::cerr << "Extending contig " << contigSoFar.seq << " by " << lastNucleotide << "\n";
#endif
            //Contig newCtg(contigSoFar);
            contigSoFar.seq += lastNucleotide;

            str_uint_map_t::const_iterator ct = kmerFreqTable_.find(s.first);
            assert(ct != kmerFreqTable_.end());
            contigSoFar.avgCoverage += ct->second;
            ++contigSoFar.numKmers;
            seenVertices[s.first] = true;
            /*std::string newSeq =*/ doPathContract_(contigSoFar,seenVertices, newKmerFreqTable, newAdjList);
            //contigSoFar.seq = newSeq;
            movedForward = true;
        }
        if (!movedForward)
        {
            std::cerr << "Dead end, adding new contig " << contigSoFar.seq << std::endl;
            newKmerFreqTable[contigSoFar.seq] = contigSoFar.avgCoverage;
        }
    }
    else
    {
        // branching point in graph, add new node here
#ifdef DEBUG_DB_GRAPH
        std::cerr << "Branching point with " << numNeighbours << " neighbours\n";
        std::cerr << "Adding new node " << contigSoFar.seq << " and recursing into neighbours\n";
#endif

        // update new adjacency list
        newKmerFreqTable[contigSoFar.seq] = contigSoFar.avgCoverage;

        for (AdjList::const_iterator ct = findIter->second.begin();
             ct != findIter->second.end(); ++ct)
        {
            if (seenVertices[ct->first]) continue;
            if (ct->second != AFFIX_TYPE::SUFFIX) continue;
            Contig newContig;
            newContig.seq = ct->first;
            //std::cerr << "New contig " << ct->first << "\n";
            ++newContig.numKmers;
            str_uint_map_t::const_iterator whIter = kmerFreqTable_.find(ct->first);
            assert(whIter != newKmerFreqTable.end());
            newContig.avgCoverage += whIter->second;
            ++newContig.numKmers;
            seenVertices[ct->first] = true;
            std::string newCtgSeq = doPathContract_(newContig,seenVertices,newKmerFreqTable,newAdjList);
            // update adjacency list
            //std::cerr << "Update adj : " << contigSoFar.seq << " => " << newCtgSeq << "\n";
            newAdjList[contigSoFar.seq].push_back(make_pair(newCtgSeq,ct->second));
        }
        //contigSoFar.clear();
    }
    return contigSoFar.seq;
}

void
DeBruijnGraph::dumpToFile(const std::string& outFilePrefix, const std::string& ref, const std::string& src) const
{
    //std::cerr << __PRETTY_FUNCTION__ <<  " outFilePrefix= " << outFilePrefix << "\n";

    std::ofstream outDotFile;
    std::ofstream outTxtFile;
    std::ofstream outJsonFile;
    //std::ofstream nodeTxtFile;

    std::stringstream sstr;
    sstr << outFilePrefix << ".k";
    sstr << wordLength_;
    //sstr << ".iter";
    //sstr << iteration;

    std::string outFileStem(sstr.str().c_str());
    std::string dotFileName  = outFileStem + ".dot";
    std::string jsonFileName = outFileStem + ".json";
    std::string txtFileName  = outFileStem + ".txt";
    //std::string nodeFileName = outFileStem + ".nodes";

    outDotFile.open(dotFileName.c_str());
    outJsonFile.open(jsonFileName.c_str());
    outTxtFile.open(txtFileName.c_str());
    //nodeTxtFile.open(nodeFileName.c_str());

    // use colors to show whether k-mers match the reference or not
    static const std::string noRefMatchColor("red");
    static const std::string sourceNodeColor("blue");
    static const std::string refMatchColor("gray");

    outDotFile << "graph {\n";
    outDotFile << "node [ style = filled ];\n";

    outJsonFile << "{\n";
    outJsonFile << "\"Nodes\" : [\n";

    str_uint_map_t aliasH;
    unsigned n(0);
     
    //std::cout << "Writing DB graph to file  " << "\n";
    //std::cout << "ref=" << ref << "\n";

    for (str_uint_map_t::const_iterator ct = kmerFreqTable_.begin(); ct!=kmerFreqTable_.end(); ++ct)
    {
        bool isRef(false); // indicates k-mer match with reference
        if (ref.find(ct->first) != std::string::npos) {
            isRef=true;
        } 
        bool isSrc(false);
        if (ct->first.find(src) != std::string::npos) {
            isSrc=true;
        } 

        aliasH[ct->first] = n++;
        //graphviz out
        //std::string newSeq = ct->first.substr(wordLength-1);
        outDotFile << aliasH[ct->first] << "[label=" << ct->first << " ";
        //outDotFile << aliasH[ct->first] << "[label=\"cov" << ct->second << "\" ";

        //json out
        outJsonFile << "{\n";
        outJsonFile << "\"id\" : " << aliasH[ct->first] << ",\n";
        outJsonFile << "\"label\" : \"" << ct->first << "\",\n";
        outJsonFile << "\"cov\" : \"" << ct->second << "\",\n";

        if (isSrc)
        {
            outDotFile << "color=" << sourceNodeColor << "]";
            outJsonFile << "\"color\" : \"" << sourceNodeColor <<  "\"" ;
        }
        else if (isRef)
        {
            outDotFile << "color=" << refMatchColor << "]";
            outJsonFile << "\"color\" : \"" << refMatchColor <<  "\"" ;
        }
        else
        {
            outDotFile << "color=" << noRefMatchColor << "]";
            outJsonFile << "\"color\" : \"" << noRefMatchColor <<  "\"\n";
        }

        if (is_last(ct,kmerFreqTable_))
        {
            outJsonFile << "\n}";
        }
        else
        {
            outJsonFile << "\n},\n";
        }
        outDotFile << "\n";
        // txt output
        outTxtFile << "VT " << aliasH[ct->first] << " " << ct->first << " " << ct->second << "\n";
    }
    outJsonFile << "\n],\n";

    outJsonFile << "\"Edges\" : [\n";

    str_uint_map_t::const_iterator last_ct1_with_value (kmerFreqTable_.begin());
    for (str_uint_map_t::const_iterator ct1 = kmerFreqTable_.begin(); ct1 !=kmerFreqTable_.end(); ++ct1)
    {
        str_vec_map_t::const_iterator ct2 = adjList_.find(ct1->first);
        //assert(ct2 != adjList.end());
        if (ct2 == adjList_.end()) continue;
        last_ct1_with_value = ct1;
    }

    for (str_uint_map_t::const_iterator ct1 = kmerFreqTable_.begin(); ct1 !=kmerFreqTable_.end(); ++ct1)
    {
        str_vec_map_t::const_iterator ct2 = adjList_.find(ct1->first);
        //assert(ct2 != adjList.end());
        if (ct2 == adjList_.end()) continue;

        //outJsonFile << " " << aliasH[ct1->first] << " : ";
        for (auto ct3 = ct2->second.begin(); ct3!=ct2->second.end(); ++ct3)
        {
            outDotFile << aliasH[ct1->first] << " -- " << aliasH[ct3->first] << ";\n";
            outJsonFile << "{ \"from\" : " <<  aliasH[ct1->first] << ", \"to\" : " << aliasH[ct3->first] << "}";

            if(ct1 == last_ct1_with_value && std::next(ct3) == ct2->second.end()) 
            {
                outJsonFile << "\n";
            }
            else
            {
                outJsonFile << ",\n";
            }

            outTxtFile << "ED " << aliasH[ct1->first] << " " << aliasH[ct3->first] << "\n";
            /*if (is_last(ct3,ct2->second))
            {
                outJsonFile << "}\n";
            }
            else
            {
                outJsonFile << "}, ";
            }*/
        }
    }
    outDotFile << "}\n";
    outJsonFile << "]\n}\n";

    /*for (auto it : aliasH) {
        nodeTxtFile << aliasH[it.first] << " " << it.second << "\n";
    }*/

    outDotFile.close();
    outTxtFile.close();
    outJsonFile.close();
    //nodeTxtFile.close();
}

void
DeBruijnGraph::setRefKmers(const std::string& ref)
{
    std::vector<std::string> words;
    generateAllKmers(ref,wordLength_,words);

    for (auto w : words)
    {
        if (kmerFreqTable_.find(w) != kmerFreqTable_.end())
        {
            kmerRefTable_[w] = true;
        }
    }
}

} // end of namespace rumovsky


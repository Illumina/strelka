
#include "StringGraph.hh"

namespace rumovsky {

    static const unsigned maxPrefixToCorrect=50;
    /// Build assembly graph from reads
    /// \return true if successful
    bool 
    StringGraph::build(const AssemblyReadInput& reads, const AssemblyReadOutput& readInfo) {

        assert(reads.size() == readInfo.size());

        const unsigned numReads(reads.size());
        for (unsigned i(0); i<numReads;++i)
        {
            if (!readInfo[i].isUsed)
            {
                std::cout << "Adding read " << reads[i].second << "\n";
                allReads_[reads[i].second].push_back(i);
            }
        }
        // build graph
        collectWords_();
        computeOverlaps_();

        // simplify
        reduceGraph_();
        return true;
    }

    void 
    StringGraph::dumpToFile(const std::string& outFilePrefix) const {
        std::ofstream ofFile;
        ofFile.open(outFilePrefix.c_str());
        ofFile << "digraph assembly {\n";
        int nodeNum(0);

        for (NodeStoreT::const_iterator i(nodes_.begin()); i!=nodes_.end(); ++i)
        {
            for (std::list<Edge>::const_iterator j(i->edges_.begin()); j!=i->edges_.end(); ++j)
            {
                ofFile << "node_" << nodeNum << " -> node_" << j->node_ << " [label=\"" << j->label_ << "-" << j->weight_ << "\"];\n";
            }
            ++nodeNum;
        }
        ofFile << "}\n";
        ofFile.close();
    }

    void 
    StringGraph::getNeighbours(const std::string& v, AdjList& adj) const {
        ReadStoreT::const_iterator ct = allReads_.find(v);
        if (ct == allReads_.begin()) {
            // vertex not found, do nothing
            return;
        }
        for (auto n : ct->second) {
            std::string nodeSeq = readSeqs_[n];
            // TODO: this is not great design, need to revisit
            adj.push_back(std::make_pair(nodeSeq,SUFFIX));
        }
    }

    bool 
    StringGraph::containsN_( const std::string& s ) { return (s.find('N') != std::string::npos); }

    bool StringGraph::compareLabelsWithNs_( const std::string& label_small, const std::string& label_large)
    {
        int numMismatches(0),numMatches(0);
        assert(label_large.size()>=label_small.size());
        for (unsigned i(0); i<label_small.size(); ++i)
        {
            if ((label_small[i]!='N')
                &&(label_large[i]!='N'))
            {
                numMismatches+=(label_small[i]!=label_large[i]);
                numMatches+=(label_small[i]==label_large[i]);
            }
        }
        return ((numMismatches==0)&&(numMatches>0));
    }

    void
    StringGraph::collectWords_()
    {
        int readNum(0);
        std::string word;

        for (ReadStoreT::const_iterator i(allReads_.begin()); i!=allReads_.end(); ++i)
        {
    #ifdef DEBUG_WORDS
            if (i->second.size()>1) {
                std::cout << "Duplicate read " << i->first << " " << i->second.size() << "\n";
            }
    #endif
            readSeqs_.push_back(i->first);

    #ifdef DEBUG_WORDS
            std::cout << "Testing read " << readSeqs_.back() << "\n";
    #endif
            for (unsigned k(0); k<=readSeqs_.back().size()-wordLength_; k++)
            {
                word.assign(readSeqs_.back().c_str()+k,wordLength_);
                if (!containsN_(word)) {
                    //std::cout << "Adding word " << word << "\n";
                    words_[word].push_back(WordLocation(readNum,k));
                } /*else {
                    std::cout << "Skipped word " << word << "\n";
                }*/
            } // ~for
            ++readNum;
        } // ~for i

    #ifdef DEBUG_WORDS
        // print list of k-mers for debugging
        std::cout << "Collected " << words_.size()  << " words: \n";
        for (WordStoreT::iterator i(words_.begin()); i!=words_.end(); ++i)
        {
            std::cout << i->first << "\n";
            for (std::vector<WordLocation>::iterator j(i->second.begin()); j!=i->second.end(); ++j)
            {
                std::cout << j->read_ << " " << j->pos_;
            }
            std::cout << "\n";
        }
    #endif
    } // ~collectWords

    ///
    /// Computes overlap length, remembers position of Ns, counts number of mismatches in overlap
    bool
    StringGraph::getOverlap_(const int readLeft,
                             const int readRight,
                             const int offset,
                             int& overlap,
                             int& diffs,
                             AllMatchToN& allMatchesToN)
    {
        int readLeftCoord(offset), readRightCoord(0);
        bool contained(false);
        std::vector<MatchToN> matchesToN;
        //  AllMatchToN allMatchesToN;
        overlap=0;
        diffs=0;

        if ((readLeftCoord==0)&&(readSeqs_[readLeft].size()==readSeqs_[readRight].size()))
        {
            // Check for two reads identical apart from Ns - otherwise they get
            // flagged as being contained in each other and both get ignored
            // TBD compare them and keep the best quality one
            contained=(readLeft>readRight);
        }
        else if (readLeftCoord<=0)
        {
            // This means readLeft is entirely contained in readRight
            readRightCoord-=readLeftCoord;
            readLeftCoord=0;
    #ifdef DEBUG_OVERLAPS
            std::cout << "CONTAINED: " << readSeqs_[readLeft] << " in " << readSeqs_[readRight] << "\n";
    #endif
            contained=true;
        }

        //  for (int i(offset);((i<reads[readLeft].size())&&(i-offset<reads[readRight].size()));i++)
        //  readLeftCoord=offset-reads[readLeft].size()+reads[readRight].size();
        //  readRightCoord=i-offset+reads[readLeft].size()-reads[readRight].size();

        for (; ( (readLeftCoord<boost::lexical_cast<int>(readSeqs_[readLeft].size()))&&(readRightCoord<boost::lexical_cast<int>(readSeqs_[readRight].size())));
                  readLeftCoord++, readRightCoord++)
        {
            //    assert(i>=0);
            //  assert(readRightCoord>=0);
            //  if (readRightCoord>=reads[readRight].size()) break;
            //    if ((reads[readLeft].c_str()[i]!='N')&&(reads[readRight].c_str()[i-offset]!='N'))
            if (readSeqs_[readLeft].c_str()[readLeftCoord]=='N')
            {
                if (readSeqs_[readRight].c_str()[readRightCoord]!='N')
                {
                    matchesToN.push_back
                    ( MatchToN( readLeft,
                                readLeftCoord,
                                readSeqs_[readRight].c_str()[readRightCoord]));
                } // ~if
            }
            else if (readSeqs_[readRight].c_str()[readRightCoord]=='N')
            {
                matchesToN.push_back
                ( MatchToN( readRight,
                            readRightCoord,
                            readSeqs_[readLeft].c_str()[readLeftCoord]));
            }
            else
            {
                // both reads are not N
                overlap++;
                //      diffs+=(reads[readLeft].c_str()[i]!=reads[readRight].c_str()[i-offset]);
    #ifdef DEBUG_OVERLAPS
                if (readSeqs_[readLeft].c_str()[readLeftCoord]
                    !=readSeqs_[readRight].c_str()[readRightCoord])
                    printf ("DIFF: %d %d\n",readLeftCoord, readRightCoord);
    #endif


                diffs+=(readSeqs_[readLeft].c_str()[readLeftCoord]
                        !=readSeqs_[readRight].c_str()[readRightCoord]);
            }
        } // ~for

        // if overlap looks good, use base calls to correct any Ns
        // TODO: number of diffs allowed should be parameter
        if (diffs<=1)
        {
            for (std::vector<MatchToN>::iterator i(matchesToN.begin());
                 i!=matchesToN.end(); ++i)
                allMatchesToN[WordLocation(i->read_,i->pos_)][i->base_]++;
            //   allMatchesToN[i->read_][i->pos_][i->base_]++;
        }
        return contained;
    } // ~getOverlap

    int
    StringGraph::countHits_(const unsigned readNum)
    {
        int numHits(0);
        std::string word;
        for (unsigned k(0); k<=readSeqs_[readNum].size()-wordLength_; k++)
        {
            word.assign(readSeqs_[readNum].c_str()+k,wordLength_);
            for (std::vector<WordLocation>::iterator j(words_[word].begin());
                 j!=words_[word].end(); ++j)
            {
                if ((j->pos_<=k)&&(j->read_!=readNum))
                {
                    numHits++;
                } // ~if
            } // ~for j
        } // ~for k
        return numHits;
    }

    void
    StringGraph::correctNs_(const AllMatchToN& allMatchesToN)
    {
        int numCorrected(0);
        //typedef map< int, map< int, map< char, int > > > AllMatchToN;
        for (AllMatchToN::const_iterator i(allMatchesToN.begin());
             i!=allMatchesToN.end(); ++i)
        {
            if (i->second.size()==1)
            {
                numCorrected++;
    #ifdef DEBUG_CORRECT
                printf ("CORRECT %d %d to %c %d\n",i->first.read_,i->first.pos_,
                        i->second.begin()->first, i->second.begin()->second);
    #endif

                assert(readSeqs_[i->first.read_][i->first.pos_]=='N');
                readSeqs_[i->first.read_][i->first.pos_]=i->second.begin()->first;
            }
        }
        printf ("corrected %d characters\n",numCorrected);
    }

    ///
    /// Computes read overlaps using starting from joint k-mers as seeds
    void
    StringGraph::computeOverlaps_()
    {
        std::string word, overlap;
        int overlapSize, numDiffs, numHits;
        bool thisReadIsContained;
        AllMatchToN allMatchesToN;
        //  readNum=0;
        for (unsigned i(0); i<readSeqs_.size(); ++i)
        {
            hits_.clear();
            numHits=0;
            nodes_.push_back( Node() );
            isContained_.push_back(false);
    #ifdef DEBUG_OVERLAPS
            std::cout << i << " " << readSeqs_[i].c_str() << "\n";
    #endif
            for (unsigned k(0); k<=readSeqs_[i].size()-wordLength_; k++)
            {
                word.assign(readSeqs_[i].substr(k,wordLength_));
                for (std::vector<WordLocation>::iterator j(words_[word].begin());
                     j!=words_[word].end(); ++j)
                {
                    //  if ((j->pos_<=k)&&(j->read_!=i))
                    //  if ((j->pos_+reads[j->read_].size()<=k+reads[i].size())&&(j->read_!=i))
                    if (    ( readSeqs_[j->read_].size()-j->pos_ >= readSeqs_[i].size()-k)
                            && (j->read_!=i))
                    {
    #ifdef DEBUG_OVERLAPS
                        printf(" %d,%d",j->read_,k-j->pos_);
    #endif
                        hits_[j->read_][k-j->pos_]++;
                        numHits++;
                    } // ~if
                } // ~for j
            } // ~for k
    #ifdef DEBUG_OVERLAPS
            printf("\n");
    #endif

            for (HitStoreT::iterator j(hits_.begin());
                 j!=hits_.end(); ++j)
            {
                for (std::map<int, int>::iterator k(j->second.begin());
                     k!=j->second.end(); ++k)
                {
    #ifdef DEBUG_OVERLAPS
                    printf("read=%d offset=%d hits=%d\n",j->first, k->first, k->second);
                    printf("%s=read\n%s=match\n",readSeqs_[i].c_str(),readSeqs_[j->first].c_str());
    #endif
                    thisReadIsContained=getOverlap_( i, j->first, k->first, overlapSize, numDiffs, allMatchesToN );
    #ifdef DEBUG_OVERLAPS
                    printf("overlap=%d errors=%d\n",overlapSize, numDiffs);
    #endif
                    //  if (k->first+k->second+WordLength-1==reads[i].size())


                    if ((numDiffs<=0)&&(overlapSize>=boost::lexical_cast<int>(overlapLength_)))
                    {
                        // TBD watch out for more than 1 overlap between same 2 reads

    #ifdef DEBUG_OVERLAPS
                        printf("HIT\n");
                        printf("read=%d offset=%d hits=%d\n",j->first, k->first, k->second);
    #endif

                        //    getOverlap( reads[i], reads[j->first], k->first, overlap);
                        if (thisReadIsContained) isContained_.back()=true;

    #ifdef DEBUG_OVERLAPS
                        if (isContained_.back()==true)
                        {
                            printf ("CONTAINED marked\n");
                        }
    #endif
    //    // TC 28.4.10 - added check that overlap extends past the first
    //    // read (which may not be true if 2nd read is quality trimmed)
    //    //      if (k->first+reads[j->first].size()>reads[i].size())
                        overlapSize=readSeqs_[i].size()-k->first;
                        if ((boost::lexical_cast<int>(readSeqs_[j->first].size())>overlapSize)&&(isContained_.back()==false))
                        {

    #ifdef DEBUG_OVERLAPS
                            std::cout << "O1: " << readSeqs_[j->first].substr(0,readSeqs_[i].size()-k->first) << std::endl
                                      << "O2: " << readSeqs_[i].substr(k->first,readSeqs_[i].size()-k->first) << std::endl
                                      << "O3: " << overlap << std::endl;
    #endif

                            assert(overlapSize<=boost::lexical_cast<int>(readSeqs_[j->first].size()));

                            overlap.assign(readSeqs_[j->first],overlapSize,
                                           readSeqs_[j->first].size()-overlapSize);
    #ifdef DEBUG_OVERLAPS
                            printf (" M %s\n", overlap.c_str()); // maximal - exact overlap
    #endif
                            //cout << "i=" << i << " " << reads[i] << endl;
                            //cout << "j=" << j->first << " " << reads[j->first] << endl;
                            assert(allReads_[readSeqs_[i]].size()*allReads_[readSeqs_[j->first]].size()!=0);
                            //return;

                            // insert edge
                            nodes_.back().edges_.push_back
                            ( Edge(j->first,
                                   overlap,
                                   overlapSize));
                            //         allReads[reads[i]]*allReads[reads[j->first]]));
                        } // ~if
                    } // ~if
                    else
                    {
    #ifdef DEBUG_OVERLAPS
                        std::cout << "read=" << j->first << " offset=" << k->first << " hits=" << k->second << "\n";
    #endif
                    } // ~else
                } // ~for k
            } // ~for j
            //    sort(nodes.back().edges_.begin(),nodes.back().edges_.end());
            nodes_.back().edges_.sort();

        } // ~for i

        correctNs_(allMatchesToN);

        for (NodeStoreT::iterator i(nodes_.begin()); i!=nodes_.end(); ++i)
        {

            for (std::list<Edge>::iterator j(i->edges_.begin()); j!=i->edges_.end(); ++j)
            {
                overlap.assign(readSeqs_[j->node_],j->weight_,
                               readSeqs_[j->node_].size()-j->weight_);

                if (overlap!=j->label_)
                {
                    std::cout << "RELABEL " << j->label_ << " to " << overlap << "\n";
                    j->label_=overlap;
                }
                //      printf ("%s %d -> node %d\n",j->label_.c_str(),j->weight_,j->node_);
                j->weight_=1;
                //      edgeNum++;
            }
            //    nodeNum++;
        }
    } // ~computeOverlaps

    void
    StringGraph::deduceMissingBases_( std::list<Edge>& edges )
    {
        std::string pref1, pref2;
        std::map<char,int> baseCounts;
        int numPrefixesFound;

    #ifdef DEBUG_DEDUCE
        for (std::list<Edge>::iterator j(edges.begin()); j!=edges.end(); ++j)
        {
            printf("DED BEF: %s %d %d\n",j->label_.c_str(),j->node_,j->weight_);
        }
    #endif

        for (unsigned i(0); i<maxPrefixToCorrect; i++)
        {
            for (std::list<Edge>::iterator j(edges.begin()); j!=edges.end(); ++j)
            {
                if (j->label_.size()<i+1) continue;
                if (j->label_.c_str()[i] == 'N')
                {
                    pref1=j->label_.substr(0,i);
    #ifdef DEBUG_DEDUCE
                    cout << "DED PREF1 " << pref1 << endl;
    #endif
                    numPrefixesFound=0;
                    baseCounts['A']=0;
                    baseCounts['C']=0;
                    baseCounts['G']=0;
                    baseCounts['T']=0;
                    for (std::list<Edge>::iterator k(edges.begin()); k!=edges.end(); ++k)
                    {
                        if (k->label_.size()<i+1) continue;
                        pref2=k->label_.substr(0,i);
                        if (pref1==pref2)
                        {
    #ifdef DEBUG_DEDUCE
                            cout << "DED PREF2 " << pref1 << endl;
    #endif
                            if (k->label_.c_str()[i]!='N')
                            {
                                baseCounts[k->label_.c_str()[i]]++;
                                numPrefixesFound++;
                            }
                        } // ~if
                    } // ~for k
                    if (numPrefixesFound>0)
                    {
    #ifdef DEBUG_DEDUCE
                        cout << "DED CNT " << baseCounts['A'] << " "<< baseCounts['C'] << " "<< baseCounts['G'] << " "<< baseCounts['T'] << " " << numPrefixesFound << endl;
    #endif
                        if (baseCounts['A']==numPrefixesFound)
                            j->label_[i]='A';
                        else if (baseCounts['C']==numPrefixesFound)
                            j->label_[i]='C';
                        else if (baseCounts['G']==numPrefixesFound)
                            j->label_[i]='G';
                        else if (baseCounts['T']==numPrefixesFound)
                            j->label_[i]='T';
    #ifdef DEBUG_DEDUCE
                        else cout << "DED no correction made" << endl;
    #endif
                    } // ~if matching prefix found
                } // ~if N found
            } // ~for j
        } // ~for i

    #ifdef DEBUG_DEDUCE
        for (std::list<Edge>::iterator j(edges.begin()); j!=edges.end(); ++j)
        {
            printf("DED AFT:  %s %d %d\n",j->label_.c_str(),j->node_,j->weight_);
        }
    #endif
    } // ~deduceMissingBases

    void
    StringGraph::reduceGraph_()
    {
        unsigned maxLabel;
        int nodeNum(0),numRemoved(0),numMarked(0),numErased(0),numContained(0);
    #ifdef DEBUG_REDUCE
        printf("reduce graph\n");
    #endif
        for (NodeStoreT::iterator i(nodes_.begin()); i!=nodes_.end(); ++i)
        {
            for (std::list<Edge>::iterator j(i->edges_.begin()); j!=i->edges_.end();)
            {
                // mark all nodes reachable from i
                if (isContained_[j->node_]==true)
                {
    #ifdef DEBUG_REDUCE
                    printf("removing link from %d to contained read %d\n",i-nodes.begin(),j->node_);
    #endif

                    numContained++;
                    j=i->edges_.erase(j);
                }
                else
                {
                    j++;
                }
            } // ~for j

        }

        for (NodeStoreT::iterator i(nodes_.begin()); i!=nodes_.end(); ++i)
        {
            deduceMissingBases_( i->edges_ );
            maxLabel=0;
    #ifdef DEBUG_REDUCE
            printf ("node %d %c %s\n", nodeNum,
                    ((isContained[nodeNum])?'Y':'N'),
                    reads[nodeNum].c_str());
    #endif
            // mark all nodes reachable from i
            for (std::list<Edge>::iterator j(i->edges_.begin()); j!=i->edges_.end(); ++j)
            {
    #ifdef DEBUG_REDUCE
                printf("EDGE %s %d %d\n",j->label_.c_str(),j->node_,j->weight_);
                //
    #endif
                //      nodes[j->node_].state_=InPlay;
                nodes_[j->node_].state_=j->weight_;
                if (j->label_.size()>maxLabel)
                    maxLabel=j->label_.size();
            }

            for (std::list<Edge>::iterator j(i->edges_.begin()); j!=i->edges_.end(); ++j)
            {
    #ifdef DEBUG_REDUCE
                printf("EDGE %s %d\n",j->label_.c_str(),j->node_);
                //
    #endif
                for (std::list<Edge>::iterator k(nodes_[j->node_].edges_.begin());
                     k!=nodes_[j->node_].edges_.end(); k++)
                {
    #ifdef DEBUG_REDUCE
                    printf("-> %s %d",k->label_.c_str(),k->node_);
                    //
    #endif

    #ifdef DEBUG_REDUCE
                    //  printf("%s %d %s %d",j->label_.c_str(),j-i->edges.begin();
                    //      k->label_.c_str(),k-i->end);
    #endif
                    if ((j->label_.size()+k->label_.size()<=maxLabel)
                        &&(nodes_[k->node_].isInPlay()))
                        //      &&(nodes[k->node_].state_==InPlay))
                    {
    #ifdef DEBUG_REDUCE
                        printf (" - REDUCING");
    #endif
                        j->weight_+=nodes_[k->node_].state_;
                        k->weight_+=nodes_[k->node_].state_;
                        nodes_[k->node_].state_=Eliminated;
                    } // ~if
    #ifdef DEBUG_REDUCE
                    printf("\n");
    #endif
                } // ~for k
            } // for j

            for (std::list<Edge>::iterator j(i->edges_.begin()); j!=i->edges_.end(); j++)
            {
                // mark all nodes reachable from i
                if (nodes_[j->node_].state_==Eliminated)
                {
                    nodes_[j->node_].state_=Vacant;
                    //  j=i->edges_.erase(j);
                    if (j->visited_==false)
                    {
                        j->visited_=true;
                        numMarked++;
                    }
                    numRemoved++;
                }
                //      else { j++; }
            } // ~for j
            nodeNum++;
        } // ~for i

        //  #ifdef DD
        for (NodeStoreT::iterator i(nodes_.begin()); i!=nodes_.end(); ++i)
        {
            for (std::list<Edge>::iterator j(i->edges_.begin()); j!=i->edges_.end();)
            {
                // mark all nodes reachable from i
                if (j->visited_==true)
                {
                    numErased++;
                    j=i->edges_.erase(j);
                }
                else
                {
                    j++;
                }
            }
        }
        //  #endif
        std::cout << "CL Removed " << numContained << " edges linking to contained reads\n" << std::endl;
        std::cout << "CL Wanted to erase " << numRemoved << " edges\n";
        std::cout << "CL Marked " << numMarked << " edges to be erased\n";
        std::cout << "CL Erased " << numErased << " edges\n";
        std::cout << "CL Total: " << numRemoved << " edges removed\n";
    } // ~reduceGraph

}

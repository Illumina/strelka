/*
 * dummy_assembler.cpp
 *
 *  Created on: Mar 12, 2015
 *      Author: mkallberg
 */

#include <assembly/dummy_assembler.hh>

dummy_assembler::dummy_assembler() {
	// TODO Auto-generated constructor stub

}

/*
 * Collected the relavant read evidence
 */
AssemblyReadInput dummy_assembler::collect_reads(int start,int end, starling_read_buffer& read_buffer){
	AssemblyReadInput reads;

	int buffer(200);
	// loop through all pile-up positions
	for (int i=start-buffer;i<end+buffer;i++){
	        read_segment_iter ri(read_buffer.get_pos_read_segment_iter(i));
	        read_segment_iter::ret_val r1;
	        read_segment_iter::ret_val r2;

	        while (true)
	        {
	            r1=ri.get_ptr();
	            if (nullptr==r1.first) break;
//
//		            // some sanity checking of read qual
//		            // TODO use command-line options for mapq cut-off (maybe shouldnt use at all?)
//		            | r1.is_invalid_realignment || !r1.is_valid()

	            if(static_cast<int>(r1.first->map_qual())>20){
	            	const read_segment r1seg(r1.first->get_segment(r1.second));
	            	const bam_seq bseq_r1(r1seg.get_bam_read());
					std::pair<int,std::string> p1(0,bseq_r1.get_string());
					reads.push_back(p1);

//						// TODO unmapped mates, where r1 indicates that the mate would be in the assembly interval
//						// need insert size distribution, likely from upstream
//						if (!nullptr==r2.first){
//							const bam_seq bseq_r2(r2.get_bam_read());
//							std::pair<int,std::string> p2(0,bseq_r2);
//							reads.push_back(p2);
//						}
	            }
				ri.next();
	        }
	}
	log_os << "Done collecting reads " << reads.size() << std::endl;
	return reads;
}

Assembly dummy_assembler::build_contigs(const AssemblyReadInput& reads, std::string& refSeq){
	Assembly contigs;
	bool verbose(false);
	const DfsAssemblerOptions opt;
	AssemblyReadOutput assembledReadInfo;
	assembledReadInfo.resize(reads.size());
	unsigned wordLength(opt.minWordLength);
	unsigned unusedReads(reads.size());

	// create graph
	rumovsky::DeBruijnGraph g(wordLength);
	static const bool addReverse(false);
	const AssemblyReadOutput readInfo;

//	readInfo.resize(reads.size());
	g.build(reads,readInfo,false);

	std::string mostFreqKmer;
	selectMostFrequentKmer(g,mostFreqKmer);

	log_os << "Setting source to most frequent k-mer " << mostFreqKmer << "\n";

	std::string src(""); // dummy source node

	// do graph pruning
	if (opt.doPruning)
		g.prune(opt.kmerPruningThreshold);

	//	log_os << "kmer freq table for DB graph class with size " <<  g.getKmerFreqTable().size() <<  "\n";

	//	log_os << "Starting greedy traversal " << "\n";
	rumovsky::GraphWalkerGreedy<rumovsky::DeBruijnGraph> walker(g);

	walker.walk(mostFreqKmer,contigs);
	//	log_os << "Graph traversal resulted in " << contigs.size() << " contigs.\n";

	const unsigned readCount(reads.size());

	// a set of read hashes; each read hash stores the starting positions of all kmers in the read
	std::vector<str_uint_map_t> readWordOffsets(readCount);

	Assembly finalContigs;
	countUnalignedReads(opt, reads, readInfo, contigs, finalContigs, readWordOffsets, unusedReads, wordLength);

	// needed for scoring eventually
	computeRead2CtgMappings(contigs,reads,readCount,wordLength,readWordOffsets,readInfo);

	log_os << "Done building contigs - found x contigs " << contigs.size() << std::endl;
	return contigs;
}

/*
 * bed_assembler.hh
 *
 *  Created on: Mar 13, 2015
 *      Author: mkallberg
 */

#ifndef C___LIB_ASSEMBLY_BED_ASSEMBLER_HH_
#define C___LIB_ASSEMBLY_BED_ASSEMBLER_HH_

#include <assembly/contigAssembler.hh>
#include "assembly/kmer/DfsAssembler.hh"
#include "assembly/kmer/DfsAssemblerOptions.hh"
#include "common/BamAddOns.hh"
#include "assembly/common/Contig.hh"

class de_bruijn_assembler: public contigAssembler {
public:
	de_bruijn_assembler();


	site_info assembleReads(int start,int end, starling_read_buffer& read_buffer,std::string& refSeq){

		// read collection and conversion to string and populate dummy read-info objects
		// TODO move to contigAssembler super class as common method for all assemblers to
		// use the same read selection method
//		int buffer(200);
//		AssemblyReadInput reads;

		// loop through all pile-up positions
//		for (int i=start-buffer;i<end+buffer;i++){
//		        read_segment_iter ri(read_buffer.get_pos_read_segment_iter(i));
//		        read_segment_iter::ret_val r1;
//		        read_segment_iter::ret_val r2;
//		        while (true)
//		        {
//		            r1=ri.get_ptr();
//		            if (nullptr==r1.first) break;
//
//		            // some sanity checking of read qual
//		            // TODO use command-line options for mapq cut-off (maybe shouldnt use at all?)
//		            if(!(static_cast<int>(r1.map_qual())<20 || r1.is_invalid_realignment || !r1.is_valid())){
//
//		            	const bam_seq bseq_r1(r1.get_bam_read());
//						std::pair<int,std::string> p1(0,bseq_r1);
//						reads.push_back(p1);
//
//						read_segment_iter::ret_val r2(r1.first->get_segment(r1.second));
//
//
//						// TODO unmapped mates, where r1 indicates that the mate would be in the assembly interval
//						// need insert size distribution, likely from upstream
//						if (!nullptr==r2.first){
//							const bam_seq bseq_r2(r2.get_bam_read());
//							std::pair<int,std::string> p2(0,bseq_r2);
//							reads.push_back(p2);
//						}
//		            }
//					ri.next();
//		        }
//		}

		// setting up some data place-holders for the assembly process
		// TO-DO the data structures can used some major cleaning here
//		Assembly as;
//		bool verbose(false);
//		DfsAssemblerOptions opt;
//		AssemblyReadOutput& assembledReadInfo;
//		Assembly& contigs;
//
//
//		assembledReadInfo.resize(reads.size());
//
//		unsigned wordLength(opt.minWordLength);
//		unsigned unusedReads(reads.size());
//
//		for (unsigned iteration(0); iteration < opt.maxAssemblyIterations; ++iteration)
//		{
////			if (unusedReads < opt.minSeedReads) break;
//
//			const unsigned lastUnusedReads(unusedReads);
//			for (; wordLength<=opt.maxWordLength; wordLength+=opt.wordStepSize)
//			{
//				if (verbose)
//				{
//					std::cerr << "Starting assembly with k=" << wordLength << " iter=" << iteration << "\n";
//					std::cerr << "opt.minWordLength=" << opt.minWordLength << " opt.maxWordLength " << opt.maxWordLength << std::endl;
//				}
//				const bool isAssemblySuccess = buildContigs(opt, reads, assembledReadInfo, wordLength,
//															contigs, unusedReads, refSeq, iteration);
//				if (verbose)
//				{
//					std::cerr << "isAssemblySuccess=" << isAssemblySuccess << std::endl;
//				}
//				if (isAssemblySuccess) break;
//			}
//
//				#ifdef DEBUG_ASBL
//						std::cerr << __PRETTY_FUNCTION__ << " iter: " << iteration << " unused readMap now: " << unusedReads << "\n";
//				#endif
//
//		        // stop if no change in number of unused reads
//		        if (unusedReads == lastUnusedReads)
//		        {
//				#ifdef DEBUG_ASBL
//							std::cerr << __PRETTY_FUNCTION__ << " Number of unused reads (" << unusedReads << ") did not change in this iteration. Stopping.\n";
//				#endif
//		        break;


		// populate site-info with assembled information /
		site_info si;
	    si.phased_ref = refSeq;
	    //const bool is_ref(max_alleles[0].first==this->reference || max_alleles[1].first==this->reference);
	    si.smod.is_block                  = false;
	    si.smod.is_unknown                = false;
//	    si.phased_alt = as[0].seq;		// for now just return the first contig
		return si;
	}
//	virtual ~bed_assembler();
};

#endif /* C___LIB_ASSEMBLY_BED_ASSEMBLER_HH_ */

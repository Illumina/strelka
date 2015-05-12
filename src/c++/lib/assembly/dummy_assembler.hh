/*
 * dummy_assembler.hh
 *
 *  Created on: Mar 12, 2015
 *      Author: mkallberg
 */

#ifndef C___LIB_ASSEMBLY_DUMMY_ASSEMBLER_HH_
#define C___LIB_ASSEMBLY_DUMMY_ASSEMBLER_HH_

#include <assembly/contigAssembler.hh>
#include "assembly/kmer/DfsAssembler.hh"
#include "assembly/kmer/DfsAssemblerOptions.hh"
#include "assembly/common/AssemblyReadInfo.hh"
#include "assembly/common/Contig.hh"
#include "common/BamAddOns.hh"
#include "assembly/common/Contig.hh"
#define DEBUG_ASSEMBLE

#ifdef DEBUG_ASSEMBLE
#include "blt_util/log.hh"
#endif




class dummy_assembler: public contigAssembler {
public:

	dummy_assembler(){};
//	virtual ~dummy_assembler();

	site_info assembleReads(int start,int end, starling_read_buffer& read_buffer,std::string& refSeq, site_info& si){
		log_os << "Im assembling " << start << "-" << end << std::endl;

		// read collection and conversion to string and populate dummy read-info objects
		// TODO move to contigAssembler super class as common method for all assemblers to
		// use the same read selection method
		int buffer(200);
		AssemblyReadInput reads;

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
//		log_os << "Done collecting " << std::endl;

		// setting up some data place-holders for the assembly process
		// TO-DO the data structures can used some major cleaning here
		Assembly as;
		bool verbose(false);
		const DfsAssemblerOptions opt;
		AssemblyReadOutput assembledReadInfo;
		Assembly contigs;
//
		assembledReadInfo.resize(reads.size());
//
		unsigned wordLength(opt.minWordLength);
		unsigned unusedReads(reads.size());

//		log_os << "Read count " << reads.size() << std::endl;
//		runDfsAssembler(opt,reads,refSeq,assembledReadInfo,contigs);

		// dummy assembly for debugging
		Contig c1;
		Contig c2;
		c1.seq = "rwegerger";
		c2.seq = "rwegerger";
		contigs.push_back(c1);
		contigs.push_back(c2);
		// populate site-info with assembled information /
	    si.phased_ref = refSeq;
	    si.smod.is_block                  = false;
	    si.smod.is_unknown                = false;
	    si.phased_alt = contigs.at(1).seq;		// for now just return the first contig
		return si;
	}

};

#endif /* C___LIB_ASSEMBLY_DUMMY_ASSEMBLER_HH_ */


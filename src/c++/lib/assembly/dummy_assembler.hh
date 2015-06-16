/*
 * dummy_assembler.hh
 *
 *  Created on: Mar 12, 2015
 *      Author: mkallberg
 */

#ifndef C___LIB_ASSEMBLY_DUMMY_ASSEMBLER_HH_
#define C___LIB_ASSEMBLY_DUMMY_ASSEMBLER_HH_

#include <assembly/contigAssembler.hh>
#include <assembly/assembly_common/AssemblyReadInfo.hh>
#include <assembly/assembly_common/BamAddOns.hh>
#include <assembly/assembly_common/Contig.hh>
#include <assembly/assembly_common/Contig.hh>
#include "assembly/kmer/DfsAssembler.hh"
#include "assembly/kmer/DfsAssemblerOptions.hh"
//#include "assembly/rumovsky_contigs.hh"
#include "assembly/graph/DeBruijnGraph.hh" // include graph
#include "assembly/kmer/GreedyAssembler.hh" // include kmer selector
#define DEBUG_ASSEMBLE

#ifdef DEBUG_ASSEMBLE
#include "blt_util/log.hh"
#endif


class dummy_assembler: public contigAssembler {
public:

	dummy_assembler();

	AssemblyReadInput collect_reads(int start,int end, starling_read_buffer& read_buffer);

	Assembly build_contigs(const AssemblyReadInput& reads, std::string& refSeq);






	site_info assembleReads(int start,int end, starling_read_buffer& read_buffer,std::string& refSeq, site_info& si){
//		log_os << "Im assembling " << start << "-" << end << std::endl;
//		log_os << "test " << std::endl;

		// read collection and conversion to string and populate dummy read-info objects
		// TODO move to contigAssembler super class as common method for all assemblers to
		// use the same read selection method

		AssemblyReadInput reads = this->collect_reads(start,end,read_buffer);
		Assembly contigs  = this->build_contigs(reads,refSeq);

		// populate site-info with assembled information /
	    si.phased_ref = contigs.at(0).seq;
	    si.n_used_calls = contigs.at(0).avgCoverage + contigs.at(1).avgCoverage;
	    si.smod.is_phased_region		  = true;
	    si.smod.is_block                  = false;
	    si.smod.is_unknown                = false;
	    si.smod.is_assembled_contig		  = true;
	    si.phased_alt = contigs.at(1).seq; // for now just return the first alternate contig

		return si;
	}

};

#endif /* C___LIB_ASSEMBLY_DUMMY_ASSEMBLER_HH_ */


/*
 * contigAssembler.hh
 *
 *  Created on: Mar 11, 2015
 *      Author: Morten Kallberg
 *      Outlines the nessacery interface for
 */

#ifndef C___LIB_ASSEMBLY_CONTIGASSEMBLER_HH_
#define C___LIB_ASSEMBLY_CONTIGASSEMBLER_HH_

#include "applications/starling/site_info_stream.hh"
#include "starling_common/starling_read_buffer.hh"

class contigAssembler{
public:
	contigAssembler(){};
//	virtual ~contigAssembler();
	virtual site_info assembleReads(int start,int end, starling_read_buffer& read_buffer,std::string& refSeq) =0;
//	{
//		site_info si;
//		return si;
//	};

	// contig assembler should own the read buffer rather than having it on assembler
	//	void register_read_buffer(starling_read_buffer& read_buffer);


};

#endif /* C___LIB_ASSEMBLY_CONTIGASSEMBLER_HH_ */

// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Starka
// Copyright (c) 2009-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

/**
 * \brief Turn assembly records into atomic calls
 *
 * \file assembly_record_processor.hh
 * \author Morten Kallberg
 * \email mkallberg@illumina.com
 *
 */


#ifndef C___LIB_ASSEMBLY_ASSEMBLY_RECORD_PROCESSOR_HH_
#define C___LIB_ASSEMBLY_ASSEMBLY_RECORD_PROCESSOR_HH_

#include "applications/starling/site_info_stream.hh"

class assembly_record_processor: public site_info_stream {
public:
    assembly_record_processor(){}

    bool add_site(site_info& si)
    {
        //TODO do something intelligent with the incoming site_info contigs starting here....
        if (si.smod.is_assembled_contig)
	{
	    si.smod.filters.set(VCF_FILTERS::AssembledRecord);
	    si.smod.is_unknown = false;
	    si.smod.is_used_covered = true;
	}
	return this->_consumer->add_site(si);
    }

    bool add_indel(const indel_info& ii)
    {
        return _consumer->add_indel(ii);
    }

    bool add_indel(const pos_t pos,
		   const indel_key ik,
		   const starling_diploid_indel_core& dindel,
		   const starling_indel_report_info& iri,
		   const starling_indel_sample_report_info& isri)
   {
     assert("unfinished method" && 0);
   }

   void flush(){}
};

#endif /* C___LIB_ASSEMBLY_ASSEMBLY_RECORD_PROCESSOR_HH_ */

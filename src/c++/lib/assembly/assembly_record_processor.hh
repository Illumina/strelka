/*
 * assembly_record_processor.hh
 *
 *  Created on: Mar 11, 2015
 *      Author: mkallberg
 */

#ifndef C___LIB_ASSEMBLY_ASSEMBLY_RECORD_PROCESSOR_HH_
#define C___LIB_ASSEMBLY_ASSEMBLY_RECORD_PROCESSOR_HH_

#include <site_info_stream.hh>

class assembly_record_processor: public site_info_stream {
public:
	assembly_record_processor(){};


	bool add_site(site_info& si){
		//TODO do something intelligent with the incoming site_info contigs starting here....
		if (si.smod.is_assembled_contig){
			si.smod.filters.set(VCF_FILTERS::AssembledRecord);
			si.smod.is_unknown = false;
			si.smod.is_used_covered = true;
		}

		return this->_consumer->add_site(si);
	}

	bool add_indel(const pos_t pos,
						  const indel_key ik,
						  const starling_diploid_indel_core& dindel,
						  const starling_indel_report_info& iri,
						  const starling_indel_sample_report_info& isri){};
	void flush(){};
};

#endif /* C___LIB_ASSEMBLY_ASSEMBLY_RECORD_PROCESSOR_HH_ */

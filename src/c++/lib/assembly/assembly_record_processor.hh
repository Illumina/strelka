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
	assembly_record_processor();
	virtual ~assembly_record_processor();
};

#endif /* C___LIB_ASSEMBLY_ASSEMBLY_RECORD_PROCESSOR_HH_ */

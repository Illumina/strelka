/*
 * site_info_stream.hh
 *
 *  Created on: Mar 10, 2015
 *      Author: Morten Kallberg
 *
 *  Abstract class specifying the site_info stream interface for
 *  linking together consumers of
 */

#include "gvcf_locus_info.hh"
//#include "starling_common/starling_base_shared.hh"

#ifndef C___SITE_INFO_STREAM_HH_
#define C___SITE_INFO_STREAM_HH_

class site_info_stream {
public:
	site_info_stream();


	virtual bool add_site(site_info& si) = 0;
	virtual bool add_indel(const pos_t pos,
						  const indel_key ik,
						  const starling_diploid_indel_core& dindel,
						  const starling_indel_report_info& iri,
						  const starling_indel_sample_report_info& isri) = 0;
	virtual void flush()  = 0;


	//TODO this needs to return a sensible value base on what is expected in starling_pos_processor around line 59
	bool is_phasing_block(){return true;}

	void clear_buffer(){
		this->_indel_buffer.clear();
		this->_site_buffer.clear();
	}

	void notify_consumer(){
		if (!this->has_listner)
			return;
		for (auto& si : this->_site_buffer)
				this->_consumer->add_site(si);
		// TODO add indel component
		//		for (auto& ii : this->_indel_buffer)
		//						this->_consumer->add_indel();
		this->clear_buffer();
	}

	void register_consumer(site_info_stream *consumer){
		this->_consumer = consumer;
		this->has_listner = true;
	}

protected:
	site_info_stream* _consumer;
	std::vector<site_info> _site_buffer;
	std::vector<indel_info> _indel_buffer;

private:
	//always make buffers available but do not populate as default
	bool has_listner = false;
};

#endif /* C___SITE_INFO_STREAM_HH_ */

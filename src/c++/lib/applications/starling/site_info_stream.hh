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
	site_info_stream(){}

	~site_info_stream(){
	  notify_consumer();
	}


	virtual bool add_site(site_info& si) = 0;

        virtual bool add_indel(const indel_info& ii) = 0;

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

        void clear_site_buffer_up_to(std::deque<site_info>::iterator & endIt){
	    _site_buffer.erase(_site_buffer.begin(),endIt);
	}

        void clear_site_buffer_to_pos(int stopPos){
	    std::deque<site_info>::iterator sit;
	    for (sit =  this->_site_buffer.begin();sit < this->_site_buffer.end();++sit)
	    {
		if ( sit->pos >= stopPos)
		{
		    break;
		}
	    }
	    clear_site_buffer_up_to(sit);
	}


        void clear_indel_buffer_up_to(std::deque<indel_info>::iterator & endIt){
	    _indel_buffer.erase(_indel_buffer.begin(),endIt);
	}

        void clear_indel_buffer_to_pos(int stopPos){
	    std::deque<indel_info>::iterator iit;
	    for (iit =  this->_indel_buffer.begin();iit < this->_indel_buffer.end();++iit)
	    {
		if ( iit->pos >= stopPos)
		{
		    break;
		}
	    }
	    clear_indel_buffer_up_to(iit);
	}


	void notify_consumer_up_to(int stopPos){
		if (!this->has_listner)
			return;

		//  for (const auto& val : _site_buffer)
	        //  {
		//    log_os << val << " ref " << val.ref << "\n";
		//  }

		std::deque<site_info>::iterator sit = _site_buffer.begin();
		std::deque<indel_info>::iterator iit = _indel_buffer.begin();
		while((sit != _site_buffer.end() && sit->pos < stopPos) || (iit != _indel_buffer.end() && iit->pos < stopPos))
		{
		    if (sit != _site_buffer.end())
		    {
                        if( iit != _indel_buffer.end())
			{
			  if ((*sit).pos >= (*iit).pos)
			  {
			      _consumer->add_indel(*iit);
			      ++iit;
			      continue;
			  }
			}
			this->_consumer->add_site(*sit);
			++sit;
			continue;
		    }
		    _consumer->add_indel(*iit);
		    ++iit;
		}
		clear_site_buffer_up_to(sit);
		clear_indel_buffer_up_to(iit);

	}


	void notify_consumer(){
	    int lastSitePos = -1;
	    if(_site_buffer.size()>0)
            {
	        lastSitePos = _site_buffer.back().pos;
	    }
	    int lastIndelPos = -1;
	    if(_indel_buffer.size()>0)
            {
	        lastIndelPos = _indel_buffer.back().pos;
	    }
	    notify_consumer_up_to(std::max(lastSitePos,lastIndelPos)+1);
	}

	void register_consumer(site_info_stream *consumer){
		this->_consumer = consumer;
		this->has_listner = true;
	}

protected:
	site_info_stream* _consumer;
	std::deque<site_info> _site_buffer;
	std::deque<indel_info> _indel_buffer;

private:
	//always make buffers available but do not populate as default
	bool has_listner = false;
};

#endif /* C___SITE_INFO_STREAM_HH_ */

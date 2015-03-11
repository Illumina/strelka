/*
 * site_info_stream.cpp
 *
 *  Created on: Mar 10, 2015
 *      Author: Morten Kallberg
 */

#include "site_info_stream.hh"

site_info_stream::site_info_stream() {
	// TODO Auto-generated constructor stub

}

site_info_stream::~site_info_stream() {
	// TODO Auto-generated destructor stub
}

void site_info_stream::clear_buffer(){
	this->_indel_buffer.clear();
	this->_site_buffer.clear();
}


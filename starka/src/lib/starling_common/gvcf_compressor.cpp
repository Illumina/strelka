/*
 * gvcfcompressor.cpp
 *
 *  Created on: Feb 21, 2014
 *  Author: Morten Kallberg
 */

#include "gvcf_compressor.hh"
#include <cassert>
#include <exception>
#include <sstream>
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <algorithm>
#include <string>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <stdlib.h>     /* atof */


#define DEBUG_GVCFCOMP

#ifdef DEBUG_GVCFCOMP
#include "blt_util/log.hh"
#endif

gvcf_compressor::gvcf_compressor() {
    this->minor_allele_loaded = false;
}

void gvcf_compressor::read_bed(const std::string& input_file, const std::string& chrom){
	    using namespace boost::algorithm;
	    std::ifstream myReadFile;
	    myReadFile.open(input_file.c_str());
	    std::string output;
	    if (myReadFile.is_open()) {
	        while (!myReadFile.eof()) {
	            std::getline (myReadFile,output);
//	            boost::replace_all(output, "chr", "");
	            std::vector<std::string> tokens;
//	            std::cout  << output << "\n";
	            split(tokens, output, is_any_of("\t")); // tokenize string
//	            //case new model
	            if (tokens.size()>3 && tokens.at(0)==chrom){
//	                std::cout  << tokens.at(0) << "\n";
//	                int my_pos = atoi( tokens.at(1).c_str() );
//	                std::cout  << pos << "\n";
//	                std::cout  << my_pos << "\n";
	                this->chr_to_pos[tokens.at(0)][atoi( tokens.at(1).c_str() )] = true;
	            }
	        }
	   }
	   this->my_chrom = chrom;
	   this->minor_allele_loaded = true;
}

bool gvcf_compressor::is_minor_allele_site(const std::string& chr, const int pos){
    chrposmap::iterator it = this->chr_to_pos.find(chr);
//    log_os << "checking for chr " << chr << " pos " << pos << "\n";

    if(it != this->chr_to_pos.end())
    {
//        log_os << "found chr " << "\n";
        posmap::iterator it2 = this->chr_to_pos[chr].find(pos);
        if(it2 != this->chr_to_pos[chr].end())
        {
//            log_os << "found pos " << "\n";
            return true;
        }
    }
    return false;
}



bool gvcf_compressor::is_site_compressable(const gvcf_options& opt, const site_info& si){

    if (! opt.is_block_compression) return false;

    if (si.dgt.is_snp) return false;

    if (si.ref!='N') {
        const double reffrac(static_cast<double>(si.known_counts[si.dgt.ref_gt]) /
                             static_cast<double>(si.n_used_calls));
        if (reffrac+opt.block_max_nonref <= 1) return false;
    }

    // check if site is in the pre-specified site that are not to be block-compressed
    if (this->minor_allele_loaded && this->is_minor_allele_site(this->my_chrom,static_cast<int>(si.pos))) return false;

    return true;
}

/*
 * gvcfcompressor.cpp
 *
 *  Created on: Feb 21, 2014
 *      Author: Morten Kallberg
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
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <stdlib.h>     /* atof */


//#define DEBUG_GVCFCOMP

#ifdef DEBUG_GVCFCOMP
#include "blt_util/log.hh"
#endif

//gvcf_compressor::gvcf_compressor() {
//	// TODO Auto-generated constructor stub
//
//}
void gvcf_compressor::read_bed(){
	using namespace boost::algorithm;
	   std::ifstream myReadFile;
	   std::cout << "reading bed \n";
	    myReadFile.open(this->minor_allele_bed.c_str());
	    std::string chr;
	    std::string pos;
	    std::string output;
	    if (myReadFile.is_open()) {
	        while (!myReadFile.eof()) {
	            std::getline (myReadFile,output);
	            std::vector<std::string> tokens;
	            std::cout  << output << "\n";
//	            split(tokens, output, is_any_of("\t")); // tokenize string
//	            //case new model
//	            if (tokens.at(0).substr(0,3)=="###") {
//	                if (pars.size()>0) {
//
//	                }
//	            }
	        }
	   }
}

// for some testing
using namespace std;
int main()
{
	gvcf_compressor comp("/home/kallberg/RefMinorAllele.bed");
	comp.read_bed();
	return(0);
}



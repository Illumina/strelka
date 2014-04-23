/*
 * gvcfcompressor.h
 *
 *  Created on: Feb 21, 2014
 *      Author: Morten Kallberg
 */

#ifndef GVCFCOMPRESSOR_H_
#define GVCFCOMPRESSOR_H_
//#include <vector>
#include <map>
#include <fstream>


class gvcf_compressor {
public:
	gvcf_compressor(
	        const std::string& bed_path) :
	        	minor_allele_bed(bed_path)
	    {}
//	virtual ~gvcf_compressor();
//    double log_odds(featuremap features, featuremap& coeffs);
    void read_bed();
//    void apply_qscore_filters(indel_info& ii, const int qscore_cut);
//    void sanity_check();
    std::string minor_allele_bed; // path to bed containing ref. is minor allele sites

};

#endif /* GVCFCOMPRESSOR_H_ */

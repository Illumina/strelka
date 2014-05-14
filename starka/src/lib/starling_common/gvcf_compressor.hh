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
#include "starling_common/gvcf_locus_info.hh"
#include "blt_common/blt_shared.hh"


typedef std::map<int, bool> posmap;
typedef std::map<std::string,posmap> chrposmap;
class gvcf_compressor {
public:
    gvcf_compressor();
    void read_bed(const std::string& input_file, const std::string& chrom);
    bool is_minor_allele_site(const std::string& chr, const int pos); //determine if the reference is the minor allele
    bool is_site_compressable(const gvcf_options& opt, const site_info& si); //determine if a site should be added to current block or if the block should be written out

//    void sanity_check();
private:
    chrposmap chr_to_pos;
    bool minor_allele_loaded;
    std::string my_chrom;
};

#endif /* GVCFCOMPRESSOR_H_ */

// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2017 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//
#include "denovo_snv_call_vcf.hh"
#include <algorithm>


void
denovo_snv_call::get_alt()
{
    /*for (unsigned i(0);i<SampleGts.size();++i){
    	std::array<unsigned,2> gt = {0,0};
    	for (unsigned chrom=0;chrom<2;chrom++){
    		if (SampleGts[i][chrom] != this->ref_gt){

    			//make sure all recorded alt alles are recorded
    			if(std::find(alts.begin(), alts.end(), SampleGts[i][chrom]) == alts.end()){
    				alts.push_back(SampleGts[i][chrom]);
    				if (alts.size()>1)
    					alt_str += ",";
    				alt_str += id_to_base(alts[alts.size()-1]);
    			}

    			//genotype sample according to ALT map
    			for (unsigned t=0;t<alts.size();t++)
    				if (alts[t]==SampleGts[i][chrom])
    					gt[chrom] = (t+1);
    		}
    	}
    	if (gt[0]>gt[1]){
    		unsigned temp = gt[1];
    		gt[1] = gt[0];
    		gt[0] = temp;
    	}
    	gts_chrom.push_back(gt);
    }

    if (alts.size()==0)
    	alt_str = ".";*/
}

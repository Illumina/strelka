// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2016 Illumina, Inc.
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

///
/// \author Chris Saunders
/// \author Morten Kallberg
///

#pragma once

#include <cstdint>
#include <array>
#include <vector>
#include <map>
#include "pedicure_vcf_locus_info.hh"
#include "blt_util/log.hh"

struct denovo_snv_call
{
    struct result_set
    {
        unsigned max_gt;
        int dsnv_qphred = 0;
        int snv_qphred = 1;
    };

    bool
    is_snv() const
    {
        return (gt_sum >0);
    }

    bool
    is_output() const
    {
        return (is_snv() || is_forced_output);
    }

    void
    consolidate_genotype(){
    	for(unsigned i=0; i<Sampleplhoods.size();i++){
    		unsigned current_min = 0;
    		unsigned sum = Sampleplhoods[i][0];
    		for (unsigned t=1; t<3; t++){
    			if (Sampleplhoods[i][t] < Sampleplhoods[i][current_min])
    				current_min = t;
    			sum += Sampleplhoods[i][t];
    		}
    		gts.push_back(current_min);
    		gqx.push_back(sum);
    		gt_sum += current_min;
    	}
    }

    void get_alt();

    unsigned ref_gt = 0;
    unsigned gt_sum = 0;
    uint8_t dsnv_tier = 0;
    bool is_forced_output = false;
    result_set rs;
    
    std::vector< std::array<float,3> > Sampleplhoods;
    std::vector< std::array<uint8_t,2> > SampleGts;
    std::vector< unsigned > gts;
    std::vector< unsigned > gqx;
    std::vector< std::array<unsigned,2> > gts_chrom;
    std::vector<uint8_t> alts;
    std::string alt_str = "";

};

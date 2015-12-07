// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Starka
// Copyright (c) 2009-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
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
    std::vector<uint8_t> alts;
    std::string alt_str = ".";

};

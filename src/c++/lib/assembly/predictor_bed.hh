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

/**
 *  \brief Predict assembly locations from bed file
 *
 *
 * \file predictor_bed.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#ifndef PREDICTOR_BED_H__
#define PREDICTOR_BED_H__

#include "predictor.hh"

#include "blt_util/RegionTracker.hh"
#include "htsapi/bed_streamer.hh"

#include <boost/algorithm/string.hpp>
#include <string>
#include <fstream>

class predictor_bed : public predictor_interface
{
public:
    predictor_bed(const char * bedfile, const char * region) : last_site (-1)
    {
        // read whole text file when not bgzipped
        if(strlen(bedfile) > 3 && strcmp(bedfile + strlen(bedfile) - 3, ".gz"))
        {
            std::ifstream in(bedfile);
            while(in.good())
            {
                std::string line;
                std::getline(in, line);
                std::vector<std::string> words;
                boost::split(words, line, boost::is_any_of("\t"));
                
                if(words.size() < 3)
                {
                    continue;
                }
                known_pos_range2 range(atoi(words[1].c_str()), atoi(words[2].c_str()));
                std::vector<std::string> contigs(words.begin()+3,words.end());
                // TODO check return code!
                this->rt.addRegion(range, contigs);
            }
        }
        else
        {
            bed_streamer *bedstr = new bed_streamer(bedfile, region);
            while ( bedstr->next() )
            {
                const bed_record *br = bedstr->get_record_ptr();
                known_pos_range2 range(br->begin, br->end);
                
                std::vector<std::string> words;
                boost::split(words, br->line, boost::is_any_of("\t"));
                
                std::vector<std::string> contigs(words.begin()+3,words.end());
                // TODO check return code!
                this->rt.addRegion(range, contigs);
            }
            delete bedstr;
        }

    }
    
    virtual void add_site(const site_info & si)
    {
        this->add_site(si.pos);
    }

    virtual void add_indel(const indel_info & ii)
    {
        int begpos = ii.pos;
        int endpos = ii.pos + ii.ref_length() - 1;
        this->add_site(std::max(begpos, endpos));
    }

    /** keep extending the current block of variants
        (this tells assembler to not forward things on to 
         the gVCF aggregator)
     */
    virtual bool keep_extending() 
    {
        return rt.isInRegion(last_site);
    }

    /**
     * @brief Return ranges to assemble
     * 
     * Last range will be [-1, -1)
     */
    std::pair<known_pos_range2, std::string> next_range()
    {
        std::pair<known_pos_range2, std::string> pr{known_pos_range2(-1, -1), ""};
        if(!positions.empty())
        {
            pr = positions.front();
            positions.pop_front();
        }
        return pr;
    }

protected:
    /** add sites / update positions */
    void add_site(int pos) 
    {
        if(last_site >= pos)
        {
            return;
        }

        int last_site_before = last_site;
        last_site = pos;
        
        boost::optional< known_pos_range2 > pl = rt.isInRegion(last_site_before);
        if(pl && !this->keep_extending())
        {
            rt.removeToPos(last_site);
            std::string p;
            boost::optional< std::vector<std::string> > pr = rt.isPayloadInRegion(last_site_before);
            if(pr)
            {
                for(std::string const & s : *pr)
                {
                    if(p.size()) { p += ","; }
                    p += s;
                }
            }
            positions.push_back(std::make_pair(*pl, p));
        }
    }

private:
    int last_site;
    std::list< std::pair<known_pos_range2, std::string> > positions;

    RegionPayloadTracker< std::vector< std::string > >  rt;
};

#endif /* PREDICTOR_BED_H__ */

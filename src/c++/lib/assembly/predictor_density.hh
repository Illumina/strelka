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
/*
 *  predictor_density.hh
 *
 *  Assembly predictor by variant density
 *
 *  Created on: Sep 10, 2013
 *  Author: Morten Kallberg / Peter Krusche / Aaron Halpern
 */

#ifndef PREDICTOR_DENSITY_H__
#define PREDICTOR_DENSITY_H__

#include "predictor.hh"

// #define DEBUG_PREDICTOR_DENSITY

class predictor_density : public predictor_interface
{
public:
    predictor_density(int _lengthToBreakRegion=20, int _varCountCutoff=3) : 
        lengthToBreakRegion(_lengthToBreakRegion),
        varCountCutoff(_varCountCutoff),
        beginNotFullyProcessed(-1),
        beginPossibleAssemblyRegion(-1),
        numVarsInPossibleAssemblyRegion(0),
        lastVarPos(-1),
        lastPos(-1) {}

    void add_site( int pos, bool isVar )
    {
        if (beginNotFullyProcessed == -1)
        {
            beginNotFullyProcessed = pos;
        }

        if (isVar)
        {
            if (beginPossibleAssemblyRegion == -1)
            {
                beginPossibleAssemblyRegion = pos;
            }
            ++numVarsInPossibleAssemblyRegion;
            if(lastVarPos < pos)
            {
                lastVarPos = pos;
            }
        }
        if (lastPos < pos) // test this in case an indel already pushed us past here
        {
            lastPos = pos;
        }

#ifdef DEBUG_PREDICTOR_DENSITY
        if ( isVar )
            std::cerr << " site addition: "
                << beginNotFullyProcessed << " : "
                <<  beginPossibleAssemblyRegion << " : "
                << numVarsInPossibleAssemblyRegion << " : "
                << lastVarPos << " : "
                << lastPos << " : "
                <<  isVar << "\n";
#endif
    }


    void add_indel( int begpos, int endpos )
    {
        if (beginNotFullyProcessed == -1)
        {
            beginNotFullyProcessed = begpos;
        }

        if (beginPossibleAssemblyRegion == -1)
        {
            beginPossibleAssemblyRegion = begpos;
        }
        ++numVarsInPossibleAssemblyRegion;
        lastVarPos = endpos;

        lastPos = endpos;

#ifdef DEBUG_PREDICTOR_DENSITY
        std::cerr << " indel addition: "
            << beginNotFullyProcessed << " : "
            <<  beginPossibleAssemblyRegion << " : "
            << numVarsInPossibleAssemblyRegion << " : "
            << lastVarPos << " : "
            << lastPos << "\n";
#endif
    }

    bool keep_extending(){
        return ( lastPos-lastVarPos < lengthToBreakRegion);
    }

    std::pair<known_pos_range2, std::string> next_range()
    {
        known_pos_range2 asmRange(-1,-1);
        if (numVarsInPossibleAssemblyRegion >= varCountCutoff)
        {


#ifdef DEBUG_PREDICTOR_DENSITY
            std::cerr << beginNotFullyProcessed << " : "
                <<  beginPossibleAssemblyRegion << " : "
                << numVarsInPossibleAssemblyRegion << " : "
                << lastVarPos << " : "
                << lastPos << "\n";
#endif

            asmRange.set_begin_pos(beginPossibleAssemblyRegion);
            asmRange.set_end_pos( lastVarPos+1 );
        }

        beginNotFullyProcessed = lastVarPos+1;
        numVarsInPossibleAssemblyRegion = 0;
        beginPossibleAssemblyRegion = -1;
        lastVarPos = -1;

        return std::make_pair(asmRange, std::string());
    }

private:
    // predictor parameters
    const int lengthToBreakRegion;           // if we see this many homref bases in a row, we know there is no assembly involving that interval
    const int varCountCutoff;                // if we see this many separate variations called without achieving lengthToBreakRegion, we should assemble

    // predictor state
    int beginNotFullyProcessed;              // positions before here have already been fully processed: assembled or not assembled, but are no longer pending
    int beginPossibleAssemblyRegion;         // bases from here on might be involved in assembly
    int numVarsInPossibleAssemblyRegion;
    int lastVarPos;                          // last position that would encourage assembly
    int lastPos;                             // last position that has been incorporated into predictor so far
};

#endif /* PREDICTOR_DENSITY_H__ */

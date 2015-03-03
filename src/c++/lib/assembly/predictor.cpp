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
 * Codon_phaser.cpp
 *
 *  Created on: Sep 10, 2013
 *  Author: Morten Kallberg
 */

#include "predictor.hh"

bool
predictor::keep_extending()
{
    //overlap with bed track
//    known_pos_range2 range(as.block_start,as.block_end);
//    if (this->rt.isInRegion(range))
//        return true;
    //more assembly condition based on reference, and other buffered metrics


    //do not assemble
    return false;
}


bool
predictor::do_assemble()
{
    //overlap with bed track
//    known_pos_range2 range(as.block_start,as.block_end);
//    if (this->rt.isInRegion(range))
//        return true;
    //more assembly condition based on reference, and other buffered metrics


    //do not assemble
    return false;
}

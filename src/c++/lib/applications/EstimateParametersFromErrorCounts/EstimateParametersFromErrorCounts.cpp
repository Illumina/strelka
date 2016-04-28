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
///

#include "EstimateParametersFromErrorCounts.hh"
#include "EPECOptions.hh"
#include "errorAnalysis/SequenceErrorCounts.hh"

#include "model1.hh"
#include "modelVariantAndIndyError.hh"
#include "modelVariantAndIndyErrorNoOverlap.hh"
#include "modelVariantAndTriggerMixError.hh"
#include "modelVariantAndBetaBinomialError.hh"

#include <iostream>



static
void
runEPEC(
    const EPECOptions& opt)
{
    SequenceErrorCounts counts;
    counts.load(opt.countsFilename.c_str());

    if (opt.modelType == MODEL_TYPE::INDEL)
    {
        if      (opt.modelIndex == 1)
        {
            model1(counts);
        }
        else if (opt.modelIndex == 2)
        {
            modelVariantAndIndyError(counts);
        }
        else if (opt.modelIndex == 3)
        {
            modelVariantAndTriggerMixError(counts);
        }
        else if (opt.modelIndex == 4)
        {
            modelVariantAndIndyErrorNoOverlap(counts);
        }
        else if (opt.modelIndex == 5)
        {
            modelVariantAndBetaBinomialError(counts);
        }
        else
        {
            std::cerr << "Unknown Indel Model\n";
        }
    }
    else
    {
        {
            std::cerr << "Unknown Model Type\n";
        }
    }
}



void
EstimateParametersFromErrorCounts::
runInternal(int argc, char* argv[]) const
{
    EPECOptions opt;

    parseEPECOptions(*this,argc,argv,opt);
    runEPEC(opt);
}

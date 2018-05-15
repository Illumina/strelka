//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2018 Illumina, Inc.
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


#include "blt_util/log.hh"
#include "EstimateVariantErrorRates.hh"
#include "EPECOptions.hh"
#include "IndelModelProduction.hh"

#include <iostream>



static
void
runEPEC(
    const EPACOptions& opt)
{
    SequenceAlleleCounts counts;
    counts.load(opt.countsFilename.c_str());


    IndelModelProduction indelModelProduction(counts, opt.thetaFilename, opt.outputFilename);
    indelModelProduction.estimateIndelErrorRates();
    if (indelModelProduction.checkEstimatedModel())
    {
        indelModelProduction.exportIndelErrorModelJson();
    }
    else
    {
        log_os << "WARNING: In EstimateVariantErrorRates on '" << counts.getSampleName()
               <<"', checkEstimatedModel() failed. Using  '" << opt.fallbackFilename << "' instead\n";
        indelModelProduction.exportModelUsingInputJson(opt.fallbackFilename);
    }
}



void
EstimateVariantErrorRates::
runInternal(int argc, char* argv[]) const
{
    EPACOptions opt;

    parseEPACOptions(*this, argc, argv, opt);
    runEPEC(opt);
}

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


#include <iostream>

#include "blt_util/log.hh"
#include "EstimateVariantErrorRates.hh"
#include "EPECOptions.hh"
#include "IndelModelProduction.hh"




static
void
runEPEC(
    const EPECOptions& opt)
{
    SequenceErrorCounts counts;
    counts.load(opt.countsFilename.c_str());


    IndelModelProduction indelModelProduction(counts, opt.thetaFilename, opt.outputFilename);
    indelModelProduction.estimateIndelErrorRates();
    if (!indelModelProduction.checkEstimatedModel())
    {
        log_os << "WARNING: In EstimateVariantErrorRates on '" << counts.getSampleName() <<"', checkEstimatedModel() failed. Using  '" << opt.fallbackFilename << "' instead\n";
        indelModelProduction.exportModelUsingInputJson(opt.fallbackFilename);
        return;
    }
    indelModelProduction.exportModel();

}



void
EstimateVariantErrorRates::
runInternal(int argc, char* argv[]) const
{
    EPECOptions opt;

    parseEPECOptions(*this,argc,argv,opt);
    runEPEC(opt);
}

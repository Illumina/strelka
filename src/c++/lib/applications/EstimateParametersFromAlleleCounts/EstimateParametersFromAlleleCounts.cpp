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

#include "EstimateParametersFromAlleleCounts.hh"
#include "EPACOptions.hh"
#include "errorAnalysis/SequenceAlleleCounts.hh"

#include "indelModel1.hh"
#include "indelModelVariantAndBetaBinomialError.hh"
#include "indelModelVariantAndIndyError.hh"
#include "indelModelVariantAndIndyErrorNoOverlap.hh"
#include "snvModel1.hh"
#include "snvModelVariantAndIndyError.hh"
#include "indelModelVariantAndBinomialMixtureError.hh"
#include "indelModelVariantAndBinomialMixtureErrorNoOverlap.hh"
#include "snvModelVariantAndBinomialMixtureError.hh"

#include <iostream>



static
void
runEPEC(
    const EPACOptions& opt)
{
    SequenceAlleleCounts counts;
    counts.load(opt.countsFilename.c_str());

    if (opt.modelType == MODEL_TYPE::INDEL)
    {
        if      (opt.modelIndex == 1)
        {
            indelModel1(counts);
        }
        else if (opt.modelIndex == 2)
        {
            indelModelVariantAndIndyError(counts);
        }
        else if (opt.modelIndex == 3)
        {
            indelModelVariantAndBinomialMixtureError(counts);
        }
        else if (opt.modelIndex == 4)
        {
            indelModelVariantAndIndyErrorNoOverlap(counts);
        }
        else if (opt.modelIndex == 5)
        {
            indelModelVariantAndBinomialMixtureErrorNoOverlap(counts);
        }
        else if (opt.modelIndex == 6)
        {
            indelModelVariantAndBetaBinomialError(counts);
        }
        else
        {
            std::cerr << "Unknown Indel Model\n";
        }
    }
    else if (opt.modelType == MODEL_TYPE::SNV)
    {
        if     (opt.modelIndex == 1)
        {
            snvModel1(counts);
        }
        else if (opt.modelIndex == 2)
        {
            snvModelVariantAndIndyError(counts);
        }
        else if (opt.modelIndex == 3)
        {
            snvModelVariantAndBinomialMixtureError(counts);
        }
        else
        {
            std::cerr << "Unknown SNV Model\n";
        }
    }
    else
    {
        std::cerr << "Unknown Model Type\n";
    }
}



void
EstimateParametersFromAlleleCounts::
runInternal(int argc, char* argv[]) const
{
    EPACOptions opt;

    parseEPACOptions(*this, argc, argv, opt);
    runEPEC(opt);
}

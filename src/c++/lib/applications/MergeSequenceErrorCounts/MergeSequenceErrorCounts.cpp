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

///
/// \author Chris Saunders
///

#include "MergeSequenceErrorCounts.hh"
#include "MSECOptions.hh"
#include "errorAnalysis/SequenceErrorCounts.hh"

#include "common/OutStream.hh"



static
void
runMSEC(const MSECOptions& opt)
{
    {
        // early test that we have permission to write to output file(s)
        OutStream outs(opt.outputFilename);
    }

    SequenceErrorCounts mergedCounts;
    bool isFirst(true);
    for (const std::string& countsFilename : opt.countsFilename)
    {
        if (isFirst)
        {
            mergedCounts.load(countsFilename.c_str());
            isFirst=false;
        }
        else
        {
            SequenceErrorCounts inputCounts;
            inputCounts.load(countsFilename.c_str());
            mergedCounts.merge(inputCounts);
        }
    }

    mergedCounts.save(opt.outputFilename.c_str());
}



void
MergeSequenceErrorCounts::
runInternal(int argc, char* argv[]) const
{
    MSECOptions opt;

    parseMSECOptions(*this,argc,argv,opt);
    runMSEC(opt);
}

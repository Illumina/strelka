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

#include "MergeSequenceAlleleCounts.hh"
#include "MSACOptions.hh"
#include "errorAnalysis/SequenceAlleleCounts.hh"

#include "common/OutStream.hh"



static
void
runMSAC(const MSACOptions& opt)
{
    {
        // early test that we have permission to write to output file(s)
        OutStream outs(opt.outputFilename);
    }

    SequenceAlleleCounts mergedCounts;
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
            SequenceAlleleCounts inputCounts;
            inputCounts.load(countsFilename.c_str());
            mergedCounts.merge(inputCounts);
        }
    }

    mergedCounts.save(opt.outputFilename.c_str());
}



void
MergeSequenceAlleleCounts::
runInternal(int argc, char* argv[]) const
{
    MSACOptions opt;

    parseMSACOptions(*this, argc, argv, opt);
    runMSAC(opt);
}

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

///
/// \author Chris Saunders
///

#include "MergeRunStats.hh"
#include "MRSOptions.hh"
#include "appstats/RunStats.hh"
#include "blt_util/log.hh"
#include "common/OutStream.hh"



static
void
runMRS(const MRSOptions& opt)
{
    {
        // early test that we have permission to write to output file(s)
        OutStream outs(opt.outputFilename);
        if (! opt.reportFilename.empty())
        {
            OutStream reps(opt.reportFilename);
        }
    }

    RunStats mergedStats;
    bool isFirst(true);
    for (const std::string& statsFilename : opt.statsFilename)
    {
        if (isFirst)
        {
            mergedStats.load(statsFilename.c_str());
            isFirst=false;
        }
        else
        {
            RunStats inputStats;
            inputStats.load(statsFilename.c_str());
            mergedStats.merge(inputStats);
        }

    }

    mergedStats.save(opt.outputFilename.c_str());
    if (! opt.reportFilename.empty())
    {
        mergedStats.report(opt.reportFilename.c_str());
    }
}



void
MergeRunStats::
runInternal(int argc, char* argv[]) const
{
    MRSOptions opt;

    parseMRSOptions(*this,argc,argv,opt);
    runMRS(opt);
}

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

#include "GetRegionDepth.hh"
#include "RegionDepthOptions.hh"
#include "ReadRegionDepthUtil.hh"

#include "blt_util/log.hh"
#include "common/OutStream.hh"

#include <cstdlib>

#include <iomanip>
#include <iostream>



static
void
getRegionDepth(const RegionDepthOptions& opt)
{
    // check that we have write permission on the output file early:
    {
        OutStream outs(opt.outputFilename);
    }

    double regionDepth(readRegionDepthFromAlignment(opt.referenceFilename, opt.alignmentFilename, opt.regions));

    OutStream outs(opt.outputFilename);
    std::ostream& os(outs.getStream());

    os << opt.alignmentFilename << "\t" << std::fixed << std::setprecision(4) << regionDepth << "\n";

}


void
GetRegionDepth::
runInternal(int argc, char* argv[]) const
{
    RegionDepthOptions opt;

    parseRegionDepthOptions(*this,argc,argv,opt);
    getRegionDepth(opt);
}

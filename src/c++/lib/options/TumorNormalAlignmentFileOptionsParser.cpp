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

#include "options/optionsUtil.hh"
#include "options/TumorNormalAlignmentFileOptionsParser.hh"

#include <set>



boost::program_options::options_description
getOptionsDescription(
    TumorNormalAlignmentFileOptions& /*opt*/)
{
    namespace po = boost::program_options;
    po::options_description desc("alignment-files");
    desc.add_options()
    ("normal-align-file", po::value<AlignmentFileOptions::files_t>(),
     "normal sample alignment file in BAM or CRAM format (exactly one file required)")
    ("tumor-align-file", po::value<AlignmentFileOptions::files_t>(),
     "tumor sample alignment file in BAM or CRAM format (exactly one file required)")
    ;
    return desc;
}


void
parseOptions(
    const boost::program_options::variables_map& vm,
    TumorNormalAlignmentFileOptions& opt)
{
    // paste together tumor and normal:
    {
        AlignmentFileOptions::files_t normal;
        AlignmentFileOptions::files_t tumor;
        if (vm.count("normal-align-file"))
        {
            normal=(boost::any_cast<AlignmentFileOptions::files_t>(vm["normal-align-file"].value()));
        }
        if (vm.count("tumor-align-file"))
        {
            tumor=(boost::any_cast<AlignmentFileOptions::files_t>(vm["tumor-align-file"].value()));
        }
        opt.alignmentFilenames = normal;
        opt.alignmentFilenames.insert(opt.alignmentFilenames.end(),
                                      tumor.begin(),
                                      tumor.end());
        opt.isAlignmentTumor.clear();
        opt.isAlignmentTumor.resize(normal.size(), false);
        opt.isAlignmentTumor.resize(opt.alignmentFilenames.size(), true);
    }
}

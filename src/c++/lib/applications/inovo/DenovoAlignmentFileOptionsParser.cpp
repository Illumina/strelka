// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

///
/// \author Chris Saunders
///

#include "DenovoAlignmentFileOptionsParser.hh"
#include "options/optionsUtil.hh"

#include <set>


typedef std::vector<std::string> files_t;



boost::program_options::options_description
getOptionsDescription(
    DenovoAlignmentFileOptions& /*opt*/)
{
    namespace po = boost::program_options;
    po::options_description desc("Sample-alignment-files");
    desc.add_options()
    ("proband-align-file", po::value<files_t>(),
     "proband alignment file in BAM or CRAM format (required)")
    ("parent-align-file", po::value<files_t>(),
     "parent alignment file in BAM or CRAM format (may be specified no more than twice)")
    ("sibling-align-file", po::value<files_t>(),
     "sibling alignment file in BAM or CRAM format (may be specified multiple times)")
    ;
    return desc;
}



static
void
addSampleGroup(
    const boost::program_options::variables_map& vm,
    const char* vmkey,
    const INOVO_SAMPLETYPE::index_t stype,
    DenovoAlignmentFileOptions& opt)
{
    if (! vm.count(vmkey)) continue;

    files_t tmpfiles=(boost::any_cast<files_t>(vm[vmkey].value()));
    for (const auto& tfile : tmpfiles)
    {
        opt.alignmentFilename.push_back(tfile);
        opt.alignmentSampleInfo.push_back(SampleInfo{stype});
    }
}



bool
parseOptions(
    const boost::program_options::variables_map& vm,
    DenovoAlignmentFileOptions& opt,
    std::string& errorMsg)
{
    using namespace INOVO_SAMPLETYPE;

    // combine all sample types into one sample list:
    addSampleGroup(vm, "proband-align-file", PROBAND, opt);
    addSampleGroup(vm, "parent-align-file", PARENT, opt);
    addSampleGroup(vm, "sibling-align-file", SIBLING, opt);

    errorMsg.clear();
    if (opt.alignmentSampleInfo.getTypeIndexList(PROBAND).size() != 1)
    {
        errorMsg="Must specify single proband alignment file";
    }
    else if (opt.alignmentSampleInfo.getTypeIndexList(PARENT).size() > 2)
    {
        errorMsg="Cannot provide more than two parent alignment files";
    }
    else
    {
        // check that alignment files exist, and names do not repeat
        std::set<std::string> nameCheck;
        for (std::string& afile : opt.alignmentFilename)
        {
            if (checkStandardizeInputFile(afile,"alignment file",errorMsg)) break;
            if (nameCheck.count(afile))
            {
                std::ostringstream oss;
                oss << "Repeated input alignment filename: " << afile << "\n";
                errorMsg = oss.str();
                break;
            }
            nameCheck.insert(afile);
        }
    }

    return (! errorMsg.empty());
}

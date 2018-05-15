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


#include "DumpSequenceAlleleCounts.hh"
#include "DSACOptions.hh"
#include "errorAnalysis/SequenceAlleleCounts.hh"

#include <iostream>
#include <iomanip>
#include <sstream>


using namespace IndelCounts;

static
void
makeExtendedIndelReportForContext(
    const IndelCounts::Context& context,
    const SingleSampleContextDataExportFormat& exportedContextData,
    std::ostream& os)
{
    static const std::string alt_sep(",");
    for (const SingleSampleContextObservationInfoExportFormat& contextObservationInfo : exportedContextData.data)
    {
        unsigned totalAltObservations = 0;
        std::ostringstream alts;
        std::ostringstream alt_counts;
        bool isFirst = true;
        for (unsigned altIndex(0); altIndex<INDEL_SIGNAL_TYPE::SIZE; ++altIndex)
        {
            if (contextObservationInfo.altObservations[altIndex] == 0) continue;

            if (! isFirst)
            {
                alts << alt_sep;
                alt_counts << alt_sep;
            }

            alts << INDEL_SIGNAL_TYPE::label(altIndex);
            alt_counts << contextObservationInfo.altObservations[altIndex];

            isFirst = false;

            totalAltObservations += contextObservationInfo.altObservations[altIndex];
        }

        // isFirst will remain true if we have seen no alts in this context
        if (isFirst)
        {
            alts << 0;
            alt_counts << 0;
        }

        // print the report
        {
            static const std::string sep("\t");

            os << std::setprecision(10);
            os << context <<  sep
               << GENOTYPE_STATUS::label(contextObservationInfo.variantStatus) << sep
               << alts.str() << sep
               << alt_counts.str() << sep
               << contextObservationInfo.refObservations << sep
               << totalAltObservations << sep
               << contextObservationInfo.contextInstanceCount << "\n"
               ;
        }
    }
}



static
void
extendedIndelReport(
    const Dataset& countsDataset,
    std::ostream& ros)
{
    ros << "context\tvariant_status\talts\talt_counts\tref_count\ttotal_alt\ttimes_observed\n";

    SingleSampleContextDataExportFormat exportedContextData;
    for (const auto& contextInfo : countsDataset)
    {
        const Context& context(contextInfo.first);
        const ContextData& contextData(contextInfo.second);

        contextData.exportData(exportedContextData);

        if (exportedContextData.data.empty()) continue;

        makeExtendedIndelReportForContext(context, exportedContextData, ros);
    }
}



static
void
runDSEC(
    const DSACOptions& opt)
{
    SequenceAlleleCounts counts;
    counts.load(opt.countsFilename.c_str());

    std::ostream& ros(std::cout);
    if (! opt.isExcludeBasecalls)
    {
        counts.getBasecallCounts().dump(ros);
    }
    if (! opt.isExcludeIndels)
    {
        Dataset& indelCounts(counts.getIndelCounts());
        if (!opt.isExtendedOutput)
        {
            indelCounts.dump(ros);
        }
        else
        {
            extendedIndelReport(indelCounts, ros);
        }
    }
}



void
DumpSequenceAlleleCounts::
runInternal(int argc, char* argv[]) const
{
    DSACOptions opt;

    parseDSACOptions(*this, argc, argv, opt);
    runDSEC(opt);
}

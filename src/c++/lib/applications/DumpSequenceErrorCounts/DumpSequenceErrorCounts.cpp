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


#include "DumpSequenceErrorCounts.hh"
#include "DSECOptions.hh"
#include "errorAnalysis/SequenceErrorCounts.hh"

#include <iostream>
#include <iomanip>
#include <sstream>


static
void
reportExtendedContext(
    const IndelErrorContext& context,
    const std::vector<ExportedIndelObservations>& observations,
    std::ostream& os)
{
    static const std::string alt_sep(",");
    for (const ExportedIndelObservations& obs : observations)
    {
        unsigned totalAltObservations = 0;
        std::ostringstream alts;
        std::ostringstream alt_counts;
        bool isFirst = true;
        for (unsigned altIndex(0); altIndex<INDEL_SIGNAL_TYPE::SIZE; ++altIndex)
        {
            if (obs.altObservations[altIndex] == 0) continue;

            if (! isFirst)
            {
                alts << alt_sep;
                alt_counts << alt_sep;
            }

            alts << INDEL_SIGNAL_TYPE::label(altIndex);
            alt_counts << obs.altObservations[altIndex];

            isFirst = false;

            totalAltObservations += obs.altObservations[altIndex];
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
               << GENOTYPE_STATUS::label(obs.variantStatus) << sep
               << alts.str() << sep
               << alt_counts.str() << sep
               << obs.refObservations << sep
               << totalAltObservations << sep
               << obs.observationCount << "\n"
               ;
        }
    }
}



static
void
extendedIndelReport(
    const IndelErrorCounts& counts,
    std::ostream& ros)
{
    ros << "context\tvariant_status\talts\talt_counts\tref_count\ttotal_alt\ttimes_observed\n";

    std::vector<ExportedIndelObservations> observations;
    for (const auto& contextInfo : counts)
    {
        const IndelErrorContext& context(contextInfo.first);
        const auto& data(contextInfo.second);

        data.exportObservations(observations);

        if (observations.empty()) continue;

        reportExtendedContext(context, observations, ros);
    }
}



static
void
runDSEC(
    const DSECOptions& opt)
{
    SequenceErrorCounts counts;
    counts.load(opt.countsFilename.c_str());

    std::ostream& ros(std::cout);
    if (! opt.isExcludeBasecalls)
    {
        counts.getBaseCounts().dump(ros);
    }
    if (! opt.isExcludeIndels)
    {
        IndelErrorCounts& indelCounts(counts.getIndelCounts());
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
DumpSequenceErrorCounts::
runInternal(int argc, char* argv[]) const
{
    DSECOptions opt;

    parseDSECOptions(*this,argc,argv,opt);
    runDSEC(opt);
}

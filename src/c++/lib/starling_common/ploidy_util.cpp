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

#include "ploidy_util.hh"
#include "blt_util/parse_util.hh"
#include "common/Exceptions.hh"
#include "htsapi/vcf_util.hh"

#include "boost/algorithm/string.hpp"

#include <sstream>



boost::optional<unsigned>
parsePloidyFromBed(const char* line)
{
    boost::optional<unsigned> result;

    if (line == nullptr) return result;

    unsigned tabcount(0);
    while (true)
    {
        if (*line=='\0' || *line=='\n') return result;
        if (*line=='\t') tabcount++;
        line++;
        if (tabcount>=4) break;
    }

    const char* s(line);
    const unsigned val = illumina::blt_util::parse_unsigned(s);
    if (s != line) result.reset(val);

    return result;
}



unsigned
parsePloidyFromBedStrict(const char* line)
{
    using namespace illumina::common;
    const auto ploidy = parsePloidyFromBed(line);
    if (! ploidy)
    {
        std::ostringstream oss;
        oss << "Can't parse ploidy (column 5) from bed record: '" << line << "'";
        BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
    }
    if (*ploidy > 2)
    {
        std::ostringstream oss;
        oss << "Parsed unsupported ploidy value (" << *ploidy << ") from bed record: '" << line << "'";
        BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
    }
    return *ploidy;
}



void
parsePloidyFromVcf(
    const unsigned expectedSampleCount,
    const char* line,
    known_pos_range2& range,
    std::vector<unsigned>& ploidy)
{
    using namespace illumina::common;
    std::ostringstream oss;

    std::vector<std::string> fields;
    boost::split(fields, line, boost::is_any_of("\t"));

    const unsigned minFieldCount(VCFID::FORMAT+expectedSampleCount);
    if (fields.size() <= minFieldCount)
    {
        oss << "Can't find expected number of fields (" << minFieldCount << ") in vcf ploidy record: '" << line << "'";
        BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
    }

    if (fields[VCFID::ALT] != "<CNV>")
    {
        oss << "Expecting ALT value of '<CNV>' in vcf ploidy record: '" << line << "'";
        BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
    }

    // set range:
    {
        range.set_begin_pos(illumina::blt_util::parse_int_str(fields[VCFID::POS]));

        static const std::string endPrefix("END=");

        std::vector<std::string> infoFields;
        boost::split(infoFields, fields[VCFID::INFO], boost::is_any_of(";"));
        bool isEndFound(false);
        for (const auto& infoKeyValuePair : infoFields)
        {
            if (boost::starts_with(infoKeyValuePair, endPrefix))
            {
                range.set_end_pos(illumina::blt_util::parse_int_rvalue(infoKeyValuePair.c_str() + endPrefix.size()));
                isEndFound=true;
                break;
            }
        }

        if (not isEndFound)
        {
            oss << "ERROR: can't find INFO/END value in vcf ploidy record: '" << line << "'\n";
            BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
        }
    }

    // set ploidy:
    ploidy.clear();
    {
        unsigned cnIndex(0);
        {
            std::vector<std::string> formatFields;
            boost::split(formatFields, fields[VCFID::FORMAT], boost::is_any_of(":"));
            bool isCNFound(false);
            for (const auto& formatTag : formatFields)
            {
                if (formatTag == "CN")
                {
                    isCNFound=true;
                    break;
                }
                cnIndex++;
            }

            if (not isCNFound)
            {
                oss << "Can't find FORMAT/CN entry in vcf ploidy record: '" << line << "'";
                BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
            }
        }

        for (unsigned sampleIndex(0); sampleIndex<expectedSampleCount; ++sampleIndex)
        {
            std::vector<std::string> sampleFields;
            boost::split(sampleFields, fields[VCFID::SAMPLE + sampleIndex], boost::is_any_of(":"));

            unsigned cn(2);
            if (sampleFields[cnIndex].size()==1)
            {
                const char c(sampleFields[cnIndex][0]);
                if (c == '1')
                {
                    cn = 1;
                }
                else if (c == '0')
                {
                    cn = 0;
                }
            }
            ploidy.push_back(cn);
        }
    }
}

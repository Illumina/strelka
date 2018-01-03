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

#include "snoise_info.hh"
#include "snoise_option_parser.hh"
#include "snoise_shared.hh"
#include "blt_util/log.hh"

#include <cstdlib>

#include <iostream>



#if 0
static
time_t
time_when_compiled()
{
    string datestr = __DATE__;
    string timestr = __TIME__;

    istringstream iss_date( datestr );
    string str_month;
    int day;
    int year;
    iss_date >> str_month >> day >> year;

    int month;
    if     ( str_month == "Jan" ) month = 1;
    else if ( str_month == "Feb" ) month = 2;
    else if ( str_month == "Mar" ) month = 3;
    else if ( str_month == "Apr" ) month = 4;
    else if ( str_month == "May" ) month = 5;
    else if ( str_month == "Jun" ) month = 6;
    else if ( str_month == "Jul" ) month = 7;
    else if ( str_month == "Aug" ) month = 8;
    else if ( str_month == "Sep" ) month = 9;
    else if ( str_month == "Oct" ) month = 10;
    else if ( str_month == "Nov" ) month = 11;
    else if ( str_month == "Dec" ) month = 12;
    else exit(-1);

    for ( string::size_type pos = timestr.find( ':' ); pos != string::npos; pos = timestr.find( ':', pos ) )
        timestr[ pos ] = ' ';
    istringstream iss_time( timestr );
    int hour, min, sec;
    iss_time >> hour >> min >> sec;

    tm t = {0};
    t.tm_mon = month-1;
    t.tm_mday = day;
    t.tm_year = year - 1900;
    t.tm_hour = hour;
    t.tm_min = min;
    t.tm_sec = sec;
    return mktime(&t);
}

static starling_info::time_t compile_time(time_when_compiled);
#endif



void
snoise_info::
usage(
    const char* xmessage) const
{
    std::ostream& os(log_os);

    static snoise_options default_opt;
    static const po::options_description visible(get_snoise_option_parser(default_opt));

    os <<
       "\n" << name() << " - strelka noise extractor\n"
       "\tversion: " << version() << "\n"
       "\n"
       "usage: " << name() << " [options]\n\n" << visible << "\n";

    if (xmessage)
    {
        os << "\n"
           << "******** COMMAND-LINE ERROR:: " << xmessage << " ********\n"
           << "\n";
        exit(EXIT_FAILURE);
    }

    exit(EXIT_SUCCESS);
}

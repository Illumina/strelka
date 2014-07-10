// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Starka
// Copyright (c) 2009-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

#include "countFastaBases.hh"

#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>


std::ostream& log_os(std::cerr);
std::ostream& report_os(std::cout);



static
void
usage(
    const char* pname)
{
    std::ostream& os(log_os);

    os <<
       "\n" << pname << " - count bases in a fasta file\n"
       "\n"
       "usage: \n" << pname << " < fasta_file\nor:\n" << pname << " fasta_file1 [[fasta_file2]...]\n"
       "\n"
       "Prints tab-delimited pair of known and total base counts, where known={ACGTacgt}.\n";

    exit(EXIT_FAILURE);
}



static const bool isBase[128] =
{
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    0,1,0,1, 0,0,0,1, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 1,0,0,0, 0,0,0,0, 0,0,0,0,
    0,1,0,1, 0,0,0,1, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 1,0,0,0, 0,0,0,0, 0,0,0,0
};



struct scan_result
{
    scan_result(const char* f)
        : fileName(f), _isSet(false) {}

    ~scan_result()
    {
        report();
    }

    void
    reset(const char* s,
          const char* e)
    {
        contigName.assign(s,e);
        knownCount=0;
        totalCount=0;
        _isSet=true;
    }

    void
    report()
    {
        if (! _isSet) return;
        report_os << fileName << '\t'
                  << contigName << '\t'
                  << knownCount << '\t'
                  << totalCount << '\n';
        _isSet=false;
    }

    const char* fileName;
    std::string contigName;
    unsigned knownCount;
    unsigned totalCount;
private:
    bool _isSet;
};



inline
const char*
next_word(const char* c,
          const bool is_end=false)
{
    while ((*c!='\0') && (is_end ^ static_cast<bool>(isspace(*c))))
    {
        c++;
    }
    return c;
}



static
bool
get_seq_counts(std::istream& ref_is,
               const char* fileName)
{
    static const unsigned buff_size(512);
    char buff[buff_size];

    bool is_header_wrap(false);
    unsigned line_no(0);
    scan_result sr(fileName);

    while (true)
    {
        bool is_wrap(false);
        ref_is.getline(buff,buff_size);
        if (! ref_is)
        {
            if     (ref_is.eof()) break;
            else if (ref_is.fail())
            {
                if (ref_is.bad())
                {
                    log_os << "ERROR: unexpected failure while attempting to read line " << (line_no+1) << "\n";
                    return false;
                }
                ref_is.clear();
            }
            is_wrap=true;
        }
        else
        {
            ++line_no;
        }

        if (is_header_wrap)
        {
            if (!is_wrap) is_header_wrap=false; // skip the remainder of wrapped header line
        }
        else if (buff[0] == '>')
        {
            if (is_wrap) is_header_wrap=true;
            sr.report();
            const char* ns(next_word(buff+1));
            if (*ns=='\0')
            {
                log_os << "ERROR: unexpected header format on line " << (line_no+is_wrap) << " : '" << buff << "'\n";
                return false;
            }
            const char* ne(next_word(ns+1,true));
            sr.reset(ns,ne);
        }
        else
        {
            if ((line_no+is_wrap)<=1)
            {
                log_os << "ERROR: missing fasta header\n";
                return false;
            }

            for (const char* b(buff); *b; ++b)
            {
                if ('\r' == *b) continue; //windows fasta files may still have '\r'
                ++sr.totalCount;
                if (isBase[static_cast<int8_t>(*b)]) ++sr.knownCount;
            }
        }
    }

    return true;
}



static
void
check_get_seq_counts(std::istream& is,
                     const char* file)
{
    if (get_seq_counts(is,file)) return;
    log_os << "ERROR: Failed to parse fasta file/stream: '" << file << "'\n";
    exit(EXIT_FAILURE);
}



void
countFastaBases::
runInternal(int argc, char* argv[]) const
{
    static const std::string helpStr[] = {"--help", "-help", "-h"};
    static const std::string* helpStrEnd = helpStr + sizeof(helpStr) / sizeof(helpStr[0]);

    if (argc==1) check_get_seq_counts(std::cin,"stdin");

    for (int i(1); i < argc; ++i)
    {
        if ((i==1) && (helpStrEnd != std::find(helpStr, helpStrEnd, argv[i]))) usage(name());
        std::ifstream ifs(argv[i]);
        if (!ifs)
        {
            log_os << "ERROR: Failed to open fasta file '" << argv[i] << "'\n";
            exit(EXIT_FAILURE);
        }
        check_get_seq_counts(ifs,argv[i]);
    }
}


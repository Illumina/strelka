// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

/// \file

/// \author Chris Saunders
///
#include "blt_util/blt_exception.hh"
#include "blt_util/export_line_parser.hh"

#include <cassert>
#include <cstdlib>
#include <cstring>

#include <iostream>
#include <sstream>



void
export_line_parser::
set_export_line(char* line) {
    static const char sep('\t');

    assert(line);

    char* p(line);

    _entry[0]=p;
    unsigned i(1);
    while ((*p != '\n') && (*p != '\0')) {
        if (*p == sep) {
            if (i == MAX_ENTRY) break;
            *p = '\0';
            _entry[i++] = p+1;
        }
        ++p;
    }

    *p = '\0';
    _is_set=true;

    if (i < SIZE) {
        _entry[i] = 0;
        std::ostringstream oss;
        oss << "ERROR:: Detected only " <<  i << " fields in export line where at least " << SIZE << " fields were expected.\n"
            << "\texport line: ";
        write_export_line(oss);
        oss << "\n";
        _is_set=false;
        throw blt_exception(oss.str().c_str());
    }
    _entry_count=i;
}



void
export_line_parser::
write_export_line(std::ostream& os) const {

    for (unsigned i(0); i<SIZE; ++i) {
        if (_entry[i]==0) break;
        if (i) os << '\t';
        os << _entry[i];
    }
}



void
export_line_parser::
str2i_die(const char* s,
          const char* type) const {
    std::ostringstream oss;
    oss << "ERROR:: export_line_parser failed to parse: '" << s << "' as type " << type << " in export line:\n\n";
    write_export_line(oss);
    oss << "\n";
    throw blt_exception(oss.str().c_str());
}



export_single_line_parser::
export_single_line_parser(const char* line) {
    assert(line);
    const unsigned line_size(strlen(line));
    assert(line_size<_line_buf_size);
    strncpy(_line_buf,line,_line_buf_size-1);
    set_export_line(_line_buf);
}

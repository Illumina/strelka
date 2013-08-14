// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/downloads/sequencing/licenses/>.
//

/// \file

/// \author Chris Saunders
///
#include "blt_util/blt_exception.hh"
#include "blt_util/export_stream_reader.hh"
#include "blt_util/log.hh"

#include <cassert>
#include <cstdlib>
#include <cstring>

#include <fstream>
#include <iostream>
#include <sstream>



bool
export_stream_reader::
next() {

    if (_isp==0) return false;

    _isp->getline(_line_buf,_line_buf_size);

    if (! *_isp) {
        if (_isp->eof()) { // normal eof:
            _is_line_set=false;
            return false;
        }

        std::ostringstream oss;
        oss << "ERROR: Unexpected read failure in export_stream_reader.next(). Current object state:\n";
        oss << "\texport_stream_label: " << name() << "\n";
        oss << "\tfailed attempting to read export_line_no: " << (_line_no+1) << "\n";
        if ((! _isp->bad()) && ((strlen(_line_buf)+1) == _line_buf_size)) {
            oss << "\texport line length possibly exceeded line buffer size of: " << (_line_buf_size-1) << "\n";
        }
        throw blt_exception(oss.str().c_str());
    }

    try {
        _is_line_set=true;
        _line_no++;
        _elp.set_export_line(_line_buf);
    } catch (const blt_exception&) {
        log_os << "ERROR:: Exception caught in export_stream_reader.next() Current object state:\n";
        report_state(log_os);
        throw;
    }
    return true;
}



void
export_stream_reader::
report_state(std::ostream& os) const {

    const export_line_parser* exlp(exline());

    os << "\texport_stream_label: " << name() << "\n";
    if (exlp) {
        os << "\texport_line_no: " << line_no() << "\n"
           << "\texport_line: ";
        exlp->write_export_line(os);
        os << "\n";
    } else {
        os << "\tno export_line currently set\n";
    }
}




export_file_reader::
export_file_reader(const char* filename)
    : _fisp(new std::ifstream(filename)) {

    if ( (! _fisp) || (! *_fisp) ) {
        log_os << "ERROR:: Can't open file: " << filename << "\n";
        exit(EXIT_FAILURE);
    }

    updateStream(*_fisp,filename);
}



export_file_reader::
~export_file_reader() { delete _fisp; }



void
export_file_reader::
rewind() {
    _fisp->clear();                   // forget we hit the end of file
    _fisp->seekg(0, std::ios::beg);   // move to the start of the file
    _line_no=0;
    _is_line_set=false;
}

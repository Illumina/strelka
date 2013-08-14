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
#include "blt_util/log.hh"
#include "blt_util/match_descriptor.hh"
#include "blt_util/parse_util.hh"
#include "blt_util/seq_util.hh"

#include <cassert>
#include <cctype>

#include <iostream>
#include <sstream>


enum {
    INDEL_BEGIN='^',
    INDEL_END='$'
};



static
void
read_length_error(const char* const md,
                  const unsigned read_length) {
    std::ostringstream oss;
    oss << "ERROR:: match descriptor: " << md << " is inconsistent with expected read_length: " << read_length << "\n";
    throw blt_exception(oss.str().c_str());
}



static
void
unknown_md_error(const char* const md,
                 const char* const mdptr) {

    std::ostringstream oss;
    oss << "ERROR:: can't parse match descriptor string: " << md << "\n"
        << "\tunexpected character: '" << *mdptr << "' at position: " << (mdptr-md+1) << "\n";
    throw blt_exception(oss.str().c_str());
}


// when IMPLICIT_SOFT_CLIP is defined, match descriptors are allowed
// to describe the match state for fewer cycles than are in the read,
// the remainder of the read is treated as if soft-clipped (i.e. as if
// it were not there):
//
#define IMPLICIT_SOFT_CLIP



void
expand_match_descriptor(const char* const read,
                        const char* const md,
                        char* const ref,
                        unsigned* const pos,
                        const unsigned read_length,
                        unsigned& ref_read_length) {

    using illumina::blt_util::parse_unsigned;

    static const char GAP('-');

    assert(NULL != read);
    assert(NULL != md);
    assert(NULL != ref);
    assert(NULL != pos);

    const char* mdptr(md);

    unsigned read_pos(0);
    unsigned map_pos(0);
    while (*mdptr) {

        if       (isdigit(*mdptr)) {
            const unsigned mlen(parse_unsigned(mdptr));
            if ((read_pos+mlen) > read_length) read_length_error(md,read_length);
            for (unsigned i(0); i<mlen; ++i) {
                ref[read_pos] = read[read_pos];
                pos[read_pos++] = ++map_pos;
            }

        } else if (is_valid_base(*mdptr)) {
            if ((read_pos+1) > read_length) read_length_error(md,read_length);
            if (*mdptr == read[read_pos]) {
                std::ostringstream oss;
                oss << "ERROR:: match descriptor indicates mismatch where none exists.\n"
                    << "\tread_pos: " << (read_pos+1) << "\n"
                    << "\tread: " << read << "\n"
                    << "\tmd: " << md << "\n";
                throw blt_exception(oss.str().c_str());
            }
            ref[read_pos] = *mdptr++;
            pos[read_pos++] = ++map_pos;

        } else if (*mdptr == INDEL_BEGIN) {
            mdptr++; // eat INDEL_BEGIN

            while (*mdptr != INDEL_END) {
                if       (isdigit(*mdptr)) {
                    const unsigned mlen(parse_unsigned(mdptr));
                    if ((read_pos+mlen) > read_length) read_length_error(md,read_length);
                    for (unsigned i(0); i<mlen; ++i) {
                        ref[read_pos] = GAP;
                        pos[read_pos++] = 0;
                    }

                } else if (is_valid_base(*mdptr)) {
                    ++mdptr;
                    ++map_pos;

                } else {
                    unknown_md_error(md,mdptr);

                }
            }

            mdptr++; // eat INDEL_END

        } else {
            unknown_md_error(md,mdptr);
        }
    }

    if (read_pos != read_length) {
#ifdef IMPLICIT_SOFT_CLIP
        if (read_pos < read_length) {
            while (read_pos<read_length) {
                ref[read_pos] = GAP;
                pos[read_pos++] = 0;
            }
        } else {
            read_length_error(md,read_length);
        }
#else
        read_length_error(md,read_length);
#endif
    }

    ref_read_length = map_pos;
    ref[read_pos] = 0;
}

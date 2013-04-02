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
#ifndef __EXPORT_LINE_PARSER_HH
#define __EXPORT_LINE_PARSER_HH

#include <boost/utility.hpp>

#include <cassert>
#include <cerrno>
#include <ciso646>
#include <climits>
#include <cstdlib>

#include <iosfwd>
#include <limits>
#include <string>


/// export_line_parser is only used directly by export_stream_reader
/// -- it modifies the ctor argument string
///
/// export_single_line_parser, below, can be used as a standalone
///
struct export_line_parser : private boost::noncopyable {

    const char* machine() const { return _entry[Machine]; }
    int run_number() const { return str2int(_entry[RunNumber]); }
    const char* run_number_str() const { return _entry[RunNumber]; }
    int lane() const { return str2int(_entry[Lane]); }
    const char* lane_str() const { return _entry[Lane]; }
    int tile() const { return str2int(_entry[Tile]); }
    const char* tile_str() const { return _entry[Tile]; }
    int x_coordinate() const { return str2int(_entry[XCoordinate]); }
    const char* x_coordinate_str() const { return _entry[XCoordinate]; }
    int y_coordinate() const { return str2int(_entry[YCoordinate]); }
    const char* y_coordinate_str() const { return _entry[YCoordinate]; }
    const char* index_string() const { return _entry[IndexString]; }
    int read_number() const { return str2int(_entry[ReadNumber]); }
    const char* read_number_str() const { return _entry[ReadNumber]; }
    const char* read () const { return _entry[Read]; }
    const char* quality() const { return _entry[QualityString]; }
    const char* match_chromosome() const { return _entry[MatchChromosome]; }
    const char* match_contig() const { return _entry[MatchContig]; }
    bool is_match_position() const { return (*_entry[MatchPosition]!='\0'); }
    long match_position() const { return str2long(_entry[MatchPosition]); }
    char match_strand() const { return _entry[MatchStrand][0]; }
    const char* match_descriptor() const { return _entry[MatchDescriptor]; }
    bool is_single_alignment_score() const { return (*_entry[SingleReadAlignmentScore]!='\0'); }
    int single_alignment_score() const { return str2int(_entry[SingleReadAlignmentScore]); }
    bool is_paired_alignment_score() const { return (*_entry[PairedReadAlignmentScore]!='\0'); }
    int paired_alignment_score() const { return str2int(_entry[PairedReadAlignmentScore]); }
    const char* partner_chromosome() const { return _entry[PartnerChromosome]; }
    const char* partner_contig() const { return _entry[PartnerContig]; }
    bool is_partner_offset() const { return (*_entry[PartnerOffset]!='\0'); }
    long partner_offset() const { return str2long(_entry[PartnerOffset]); }
    bool is_partner_strand() const { return (partner_strand()!='\0'); }
    char partner_strand() const { return (*_entry[PartnerStrand]); }
    bool is_passed_filter() const { return (_entry[Filtering][0] == 'Y'); }

    /// write the full export line w/o the newline character:
    void write_export_line(std::ostream& os) const;

    unsigned extra_field_count() const { return _entry_count-SIZE; }
    const char* extra_field(const unsigned i) const { assert((SIZE+i)<_entry_count); return _entry[SIZE+i]; }

protected:
    friend struct export_stream_reader;

    export_line_parser() : _is_set(false) {
        _entry[0]=0;
    }

    void
    set_export_line(char* line);

private:

    int
    str2long(const char* s,
             const char* type = "long") const {

        static const int base(10);

        errno = 0;
        char* end_ptr(0);
        const long i(strtol(s,&end_ptr,base));
        if ((errno == ERANGE && (i == LONG_MAX || i == LONG_MIN)) ||
            (errno != 0 && i == 0) ||
            (end_ptr == s)) {
            str2i_die(s,type);
        }

        return i;
    }

    int
    str2int(const char* s) const {
        const long i(str2long(s,"int"));

        if((i > std::numeric_limits<int>::max()) ||
           (i < std::numeric_limits<int>::min())) {
            str2i_die(s,"int");
        }
        return static_cast<int>(i);
    }

    void
    str2i_die(const char* s,
              const char* type) const;

    enum index_t {
        Machine,
        RunNumber,
        Lane,
        Tile,
        XCoordinate,
        YCoordinate,
        IndexString,
        ReadNumber,
        Read,
        QualityString,
        MatchChromosome,
        MatchContig,
        MatchPosition,
        MatchStrand,
        MatchDescriptor,
        SingleReadAlignmentScore,
        PairedReadAlignmentScore,
        PartnerChromosome,
        PartnerContig,
        PartnerOffset,
        PartnerStrand,
        Filtering,
        SIZE
    };

    enum { MAX_ENTRY = 50 };

    const char* _entry[MAX_ENTRY];
    unsigned _entry_count;
    bool _is_set;
};


/// this version can be used as a standalone, it does not modify the
/// ctor argument
///
struct export_single_line_parser : public export_line_parser {

    export_single_line_parser(const char* line);

private:
    enum { _line_buf_size = 2000 };
    char _line_buf[_line_buf_size];
};


#endif

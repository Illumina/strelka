// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2012 Illumina, Inc.
//
// This software is covered by the "Illumina Genome Analyzer Software
// License Agreement" and the "Illumina Source Code License Agreement",
// and certain third party copyright/licenses, and any user of this
// source file is bound by the terms therein (see accompanying files
// Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
// Illumina_Source_Code_License_Agreement.pdf and third party
// copyright/license notices).
//
//

/// \file
///
/// an efficient (and slightly unsafe) class for basic tab-delimited files, etc...
///

/// \author Chris Saunders
///
#ifndef ISTREAM_LINE_SPLITTER_HH__
#define ISTREAM_LINE_SPLITTER_HH__

#include <iosfwd>


struct istream_line_splitter {

    istream_line_splitter(std::istream& is,
                          const unsigned line_buf_size=4096,
                          const char word_seperator='\t',
                          const unsigned max_word=0)
        : _is(is)
        , _line_no(0)
        , _n_word(0)
        , _buf_size(line_buf_size)
        , _sep(word_seperator)
        , _max_word(max_word)
        , _buf(new char[_buf_size]) {

        if((0==_max_word) || (MAX_WORD_COUNT < _max_word)){
            _max_word=MAX_WORD_COUNT;
        }
    }

    ~istream_line_splitter() { if(NULL!=_buf) { delete [] _buf; _buf=NULL;} }

    unsigned
    n_word() const { return _n_word; }

    void
    dump(std::ostream& os) const;

    /// returns false for regular end of input:
    bool
    parse_line();


    enum { MAX_WORD_COUNT = 50 };
    char* word[MAX_WORD_COUNT];
private:
    std::istream& _is;
    unsigned _line_no;
    unsigned _n_word;
    unsigned _buf_size;
    char _sep;
    unsigned _max_word;
    char* _buf;
};



#if 0
{ //usage example:
    istream_line_splitter dparse(data_is);

    while(dparse.parse_line()) {
        static const unsigned col_count(46);
        if(dparse.n_word()!=col_count){
            std::ostringstream oss;
            oss << "ERROR: unexpected number of columns in paired export line:\n\n";
            dparse.dump(oss);
            throw blt_exception(oss.str().c_str());
        }

        for(unsigned i(1);(i+1)<col_count;++i){
            dparse.word[i][strlen(dparse.word[i])] = sep;
        }
        const char* nocompress_segment(dparse.word[0]);
        const char* compress_segment(dparse.word[1]);

        /// ....etc
    }
}
#endif


#endif

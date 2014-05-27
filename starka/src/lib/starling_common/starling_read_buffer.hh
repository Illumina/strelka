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
///
/// \author Chris Saunders
///

#ifndef __STARLING_READ_BUFFER_HH
#define __STARLING_READ_BUFFER_HH

#include "starling_common/starling_read.hh"

#include "boost/utility.hpp"

#include <map>
#include <set>


// Simple id incrementer, by default starling read buffer uses this
// object to increment read ids within the buffer. A reference to this
// object can also be passed to srb ctor such that read-id's will be
// unique against one or more buffers (for ie. other samples) which
// have been passed the same counter object. Obviously a single thread
// solution.
//
struct read_id_counter {
    read_id_counter() : _head_read_id(0) {}

    // provide the current id and increment:
    align_id_t next() { return _head_read_id++; }

private:
    align_id_t _head_read_id; // tracks head read_id number for next read added
};



struct read_segment_iter;


//
// Must be able to look up reads by
// (1) first alignment pos
// (2) key or
// (3) read_id_no
// (4) contig_id_no(s)
//
// multiple reads may be associated with (1) and (4), but (2) and (3)
// can produce at most a single result.
//
struct starling_read_buffer : private boost::noncopyable {

    starling_read_buffer(read_id_counter* ricp = NULL)
        : _ricp( (NULL==ricp) ? &_ric : ricp ) {}

    ~starling_read_buffer();

    /// \return <was the read added?, what is the id in the read buffer? >
    ///
    // note pos_processor is responsible for checking that the
    // position of the read is not too low -- the read_buffer itself
    // is agnostic to the data management process:
    //
    // TODO: return boost::optional
    //
    std::pair<bool,align_id_t>
    add_read_alignment(const starling_options& opt,
                       const bam_record& br,
                       const alignment& al,
                       const MAPLEVEL::index_t maplev,
                       const READ_ALIGN::index_t rat,
                       const align_id_t contig_id);

#if 1
    // adjust read segment's buffer position to new_buffer_pos,
    // and change buffer pos:
    void
    rebuffer_read_segment(const align_id_t read_id,
                          const seg_id_t seg_id,
                          const pos_t new_buffer_pos);
#endif

    read_segment_iter
    get_pos_read_segment_iter(const pos_t pos);

    // returns NULL if read_id isn't present:
    starling_read*
    get_read(const align_id_t read_id) {
        const read_data_t::iterator k(_read_data.find(read_id));
        if (k == _read_data.end()) return NULL;
        return (k->second);
    }

    // returns NULL if read_id isn't present:
    const starling_read*
    get_read(const align_id_t read_id) const {
        const read_data_t::const_iterator k(_read_data.find(read_id));
        if (k == _read_data.end()) return NULL;
        return (k->second);
    }

    void
    clear_pos(const bool is_ignore_read_names,
              const pos_t pos);

    void
    dump_pos(const pos_t pos, std::ostream& os) const;

    bool
    empty() const { return _pos_group.empty(); }

private:
    align_id_t
    next_id() const { return _ricp->next(); }

    friend struct read_segment_iter;

    //
    typedef read_key read_key_t;
    typedef std::map<align_id_t,starling_read*> read_data_t;
    typedef std::map<read_key_t, align_id_t> read_key_lup_t;
    typedef std::set<align_id_t> read_group_t;
    typedef std::pair<align_id_t,seg_id_t> segment_t;
    typedef std::set<segment_t> segment_group_t;
    typedef std::map<pos_t,segment_group_t> pos_group_t;
    typedef std::map<align_id_t,read_group_t> align_id_group_t;

    static const segment_group_t _empty_segment_group;

    // used to produce unique, sequential read ids accross multiple
    // starling_read_buffer objects (where there would typically be one
    // buffer per sample):
    //
    read_id_counter _ric; // only used if a counter isn't specified on the cmdline
    read_id_counter* _ricp;

    // read id to read data structure pointer map:
    read_data_t _read_data;

    // read name (eg. QNAME) to read id map:
    read_key_lup_t _read_key;

    // storage position to read segment id map
    //
    // note that storage position starts out as the starting position
    // of the read, however the read may be realigned without changing
    // the storage position:
    //
    pos_group_t _pos_group;

    // map from contig id to all corresponding contig reads
    align_id_group_t _contig_group;
};



// not a real iterator
//
struct read_segment_iter {

    typedef std::pair<starling_read*,seg_id_t> ret_val;

    // returns first=NULL if no read segments left:
    //
    ret_val get_ptr();

    // returns false if no more reads
    //
    bool next() {
        if (_head!=_end) _head++;
        return (_head!=_end);
    }

private:
    friend struct starling_read_buffer;
    typedef starling_read_buffer::segment_group_t::const_iterator piter;

    read_segment_iter(starling_read_buffer& buff,
                      const piter begin,
                      const piter end)
        : _buff(buff), _head(begin), _end(end) {}

    starling_read_buffer& _buff;
    piter _head;
    const piter _end;
};


#endif

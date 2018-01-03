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

/// \file
///
/// \author Chris Saunders
///

#pragma once

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
struct read_id_counter
{
    read_id_counter() : _head_read_id(0) {}

    /// provide the current id and increment:
    align_id_t next()
    {
        return _head_read_id++;
    }

private:
    align_id_t _head_read_id; // tracks head read_id number for next read added
};



struct read_segment_iter;


/// Buffer to store read information
///
/// this can be used to look up reads by
/// (1) first alignment pos
/// (2) read_id (unique id created internally)
///
/// multiple reads may be associated with (1), but (2)
/// should produce at most a single result.
///
struct starling_read_buffer : private boost::noncopyable
{
    explicit
    starling_read_buffer(read_id_counter* ricp = nullptr)
        : _ricp( (nullptr==ricp) ? &_ric : ricp ) {}

    ~starling_read_buffer();

    /// insert new read into read buffer
    ///
    /// \return what is the read's internal id in the buffer?
    ///
    // note pos_processor is responsible for checking that the
    // position of the read is not too low -- the read_buffer itself
    // is agnostic to the data management process:
    //
    align_id_t
    add_read_alignment(
        const bam_record& br,
        const alignment& inputAlignment,
        const MAPLEVEL::index_t maplev);

    /// adjust read segment's buffer position to new_buffer_pos,
    /// and change buffer pos:
    void
    rebuffer_read_segment(const align_id_t read_id,
                          const seg_id_t seg_id,
                          const pos_t new_buffer_pos);

    read_segment_iter
    get_pos_read_segment_iter(const pos_t pos);

    /// \return pointer to read, or nullptr if read_id isn't in buffer
    starling_read*
    get_read(const align_id_t read_id)
    {
        const read_data_t::iterator k(_read_data.find(read_id));
        if (k == _read_data.end()) return nullptr;
        return (k->second);
    }

    /// \return pointer to read, or nullptr if read_id isn't in buffer
    const starling_read*
    get_read(const align_id_t read_id) const
    {
        const read_data_t::const_iterator k(_read_data.find(read_id));
        if (k == _read_data.end()) return nullptr;
        return (k->second);
    }

    /// clear contents of read buffer up to and including position pos
    void
    clear_to_pos(
        const pos_t pos)
    {
        pos_group_t::iterator iter(_pos_group.begin());
        const pos_group_t::iterator end(_pos_group.upper_bound(pos));

        while (iter != end)
        {
            clear_iter(iter++);
        }
    }

    void
    dump_pos(const pos_t pos, std::ostream& os) const;

    /// clear contents of read buffer
    void
    clear()
    {
        pos_group_t::iterator iter(_pos_group.begin());
        const pos_group_t::iterator end(_pos_group.end());

        while (iter != end)
        {
            clear_iter(iter++);
        }
    }

    /// \return total reads buffered
    unsigned
    size() const
    {
        return _read_data.size();
    }

    bool
    empty() const
    {
        return _pos_group.empty();
    }

private:
    typedef std::map<align_id_t,starling_read*> read_data_t;
    typedef std::pair<align_id_t,seg_id_t> segment_t;
    typedef std::set<segment_t> segment_group_t;
    typedef std::map<pos_t,segment_group_t> pos_group_t;


    /// the input iterator will be invalidated by this method!
    void
    clear_iter(
        const pos_group_t::iterator i);

    /// \brief Generates a unique read id for each strelka process
    ///
    /// Currently this just means that we increment the read id as each read is observed when merging reads from
    /// multiple samples. Our scheme currently assumes a single thread so there is no sychronization required here.
    align_id_t
    getNextReadIndex() const
    {
        return _ricp->next();
    }

    friend struct read_segment_iter;

    static const segment_group_t _empty_segment_group;

    // used to produce unique, sequential read ids across multiple
    // starling_read_buffer objects (where there would typically be one
    // buffer per sample):
    //
    read_id_counter _ric; // only used if a counter isn't specified on the cmdline
    read_id_counter* _ricp;

    // read id to read data structure pointer map:
    read_data_t _read_data;

    // storage position to read segment id map
    //
    // note that storage position starts out as the starting position
    // of the read, however the read may be realigned without changing
    // the storage position:
    //
    pos_group_t _pos_group;
};


/// object which allows iteration over read segments
///
/// note this is not an STL iterator
struct read_segment_iter
{
    typedef std::pair<starling_read*,seg_id_t> ret_val;

    /// \return information on the read segment that the iterator current points to
    ///         ret_val.first is set to nullptr if no read segments are left
    ret_val get_ptr();

    /// advance iterator to next read segment
    /// \return false if no more reads
    bool next()
    {
        if (_head!=_end) ++_head;
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

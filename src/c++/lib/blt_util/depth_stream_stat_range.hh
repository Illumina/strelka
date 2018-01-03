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

#pragma once

#include "depth_stream_stat.hh"



/// \brief Simple on-line statistics for unsigned values
///
/// This version of stream stat requires that you specify the
/// total number of observations beforehand. All observations not
/// provided as updates() are assumed to be zero.
///
struct depth_stream_stat_range
{
    depth_stream_stat_range(const unsigned total_obs)
        : _kprior(total_obs)
        , _is_pstat(false) {}

    void update (const unsigned d)
    {
        if (_is_pstat) _is_pstat=false;
        _ustat.update(d);
    }

    unsigned sample_size() const
    {
        return get_pstat().sample_size();
    }
    unsigned nonzero() const
    {
        return get_pstat().nonzero();
    }
    double min() const
    {
        return get_pstat().min();
    }
    double max() const
    {
        return get_pstat().max();
    }
    double mean() const
    {
        return get_pstat().mean();
    }
    double variance() const
    {
        return get_pstat().variance();
    }
    double sd() const
    {
        return get_pstat().sd();
    }
    double stderror() const
    {
        return get_pstat().stderror();
    }

    friend
    std::ostream&
    operator<<(std::ostream&,const depth_stream_stat_range&);

private:

    const depth_stream_stat&
    get_pstat() const
    {
        if (! _is_pstat)
        {
            _pstat=_ustat;
            for (unsigned i(_ustat.sample_size()); i<_kprior; ++i) _pstat.update(0);
            _is_pstat=true;
        }
        return _pstat;
    }

    depth_stream_stat _ustat;
    mutable depth_stream_stat _pstat;
    unsigned _kprior;
    mutable bool _is_pstat;
};


inline
std::ostream&
operator<<(std::ostream& os,const depth_stream_stat_range& ss)
{
    return os << ss.get_pstat();
}

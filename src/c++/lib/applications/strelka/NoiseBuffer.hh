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

#include "SiteNoise.hh"
#include "blt_util/blt_types.hh"
#include "blt_util/RangeMap.hh"


struct NoiseBuffer
{
    void
    insertSiteNoise(
        const pos_t pos,
        const SiteNoise& sn)
    {
        _ndata.getRef(pos) = sn;
    }

    /// \returns nullptr for empty pos
    ///
    const SiteNoise*
    getPos(const pos_t pos) const
    {
        if (! _ndata.isKeyPresent(pos)) return nullptr;
        return &_ndata.getConstRef(pos);
    }

    void
    clear_pos(const pos_t pos)
    {
        if (_ndata.isKeyPresent(pos)) _ndata.erase(pos);
    }

    void
    clear()
    {
        _ndata.clear();
    }

    bool
    empty() const
    {
        return _ndata.empty();
    }

//    void
//    dump(std::ostream& os) const;

private:
    typedef RangeMap<pos_t,SiteNoise,ClearT<SiteNoise>> ndata_t;

    ndata_t _ndata;
};


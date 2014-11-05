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

///
/// \author Chris Saunders
///

#pragma once

#include "../../blt_util/RangeMap.hh"
#include "SiteNoise.hh"
#include "blt_util/blt_types.hh"


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


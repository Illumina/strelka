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

#include "starling_common/prog_info_base.hh"


struct snoise_info : public prog_info_base
{
    static
    const prog_info& get()
    {
        static const snoise_info vci;
        return vci;
    }

private:
    const char* name() const override
    {
        static const char NAME[] = "StrelkaNoiseExtractor";
        return NAME;
    }

    void usage(const char* xmessage = nullptr) const;

    void doc() const {}

    snoise_info() {}
    virtual ~snoise_info() {}
};

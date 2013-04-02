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
///
/// \author Chris Saunders
///

#ifndef __STARLING_INFO_HH
#define __STARLING_INFO_HH

#include "blt_util/prog_info.hh"


struct starling_info : public prog_info {

    static
    const prog_info& get() {
        static const starling_info vci;
        return vci;
    }

private:
    const char* name() const {
        static const char NAME[] = "starling";
        return NAME;
    }

    const char* version() const {
        static const char VERSION[] = "${VERSION}";
        return VERSION;
    }

    void usage(const char* xmessage = 0) const;

    void doc() const {};

    starling_info() {}
    virtual ~starling_info() {}
};

#endif

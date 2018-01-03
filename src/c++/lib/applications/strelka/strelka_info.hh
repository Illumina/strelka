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

#pragma once

#include "starling_common/prog_info_base.hh"


struct strelka_info : public prog_info_base
{
    static
    const prog_info& get()
    {
        static const strelka_info vci;
        return vci;
    }

private:
    const char* name() const override
    {
        static const char NAME[] = "strelka";
        return NAME;
    }

    void usage(const char* xmessage = 0) const override;

    void doc() const override;

    strelka_info() {}
    ~strelka_info() override {}
};

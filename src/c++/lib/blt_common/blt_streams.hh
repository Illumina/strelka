// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2016 Illumina, Inc.
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

#include "blt_common/blt_shared.hh"
#include "blt_util/prog_info.hh"
#include "htsapi/bam_util.hh"

#include <iosfwd>
#include <memory>
#include <string>


struct blt_streams
{
    blt_streams(const blt_options& opt,
                const prog_info& pinfo);

    std::ostream* report_osptr() const
    {
        return _report_osptr.get();
    }
    std::ostream* counts_osptr() const
    {
        return _counts_osptr.get();
    }
    std::ostream* nonref_test_osptr() const
    {
        return _nonref_test_osptr.get();
    }
    std::ostream* nonref_sites_osptr() const
    {
        return _nonref_sites_osptr.get();
    }

protected:

    static
    void
    open_ofstream(const prog_info& pinfo,
                  const std::string& filename,
                  const char* label,
                  std::ofstream& fos);

    /// write first few meta-data lines for a legacy blt/starling
    /// variant output file:
    ///
    static
    void
    write_file_audit(const blt_options& opt,
                     const prog_info& pinfo,
                     const char* const cmdline,
                     std::ostream& os);

    /// write the first few meta-data lines for a vcf file:
    ///
    static
    void
    write_vcf_audit(const blt_options& opt,
                    const prog_info& pinfo,
                    const char* const cmdline,
                    const bam_hdr_t& header,
                    std::ostream& os);

private:
    std::unique_ptr<std::ostream> _report_osptr;
    std::unique_ptr<std::ostream> _counts_osptr;
    std::unique_ptr<std::ostream> _nonref_test_osptr;
    std::unique_ptr<std::ostream> _nonref_sites_osptr;
};

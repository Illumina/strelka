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

#include "blt_common/blt_shared.hh"
#include "blt_util/prog_info.hh"
#include "htsapi/bam_util.hh"

#include <iosfwd>
#include <memory>
#include <string>


struct blt_streams
{
    blt_streams(const blt_options& client_opt,
                const prog_info& pinfo,
                const bool is_include_seq_name=false);

    std::ostream* report_osptr() const
    {
        return _report_osptr.get();
    }
    std::ostream* counts_osptr() const
    {
        return _counts_osptr.get();
    }
    std::ostream* bsnp_diploid_osptr() const
    {
        return _bsnp_diploid_osptr.get();
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
                  const bool is_clobber,
                  std::ofstream& fos);

    // write first few meta-data lines for a legacy blt/starling
    // variant output file:
    //
    static
    void
    write_file_audit(const blt_options& opt,
                     const prog_info& pinfo,
                     const char* const cmdline,
                     std::ostream& os);

    // write the first few meta-data lines for a vcf file:
    //
    static
    void
    write_vcf_audit(const blt_options& opt,
                    const prog_info& pinfo,
                    const char* const cmdline,
                    const bam_header_t* const header,
                    std::ostream& os);

private:
    std::unique_ptr<std::ostream> _report_osptr;
    std::unique_ptr<std::ostream> _counts_osptr;
    std::unique_ptr<std::ostream> _bsnp_diploid_osptr;
    std::unique_ptr<std::ostream> _nonref_test_osptr;
    std::unique_ptr<std::ostream> _nonref_sites_osptr;
};

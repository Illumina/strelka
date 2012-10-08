// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2012 Illumina, Inc.
//
// This software is covered by the "Illumina Genome Analyzer Software
// License Agreement" and the "Illumina Source Code License Agreement",
// and certain third party copyright/licenses, and any user of this
// source file is bound by the terms therein (see accompanying files
// Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
// Illumina_Source_Code_License_Agreement.pdf and third party
// copyright/license notices).
//
//

/// \file

/// \author Chris Saunders
///
#ifndef __BLT_STREAMS_HH
#define __BLT_STREAMS_HH

#include "blt_common/blt_shared.hh"
#include "blt_util/prog_info.hh"

#include <iosfwd>
#include <memory>
#include <string>


struct blt_streams {

    blt_streams(const blt_options& client_opt,
                const prog_info& pinfo,
                const bool is_include_seq_name=false);

    ~blt_streams();

    std::ostream& report_os() const { return _report_os; }
    std::ostream* counts_osptr() const { return _counts_osptr.get(); }
    std::ostream* bsnp_diploid_osptr() const { return _bsnp_diploid_osptr.get(); }
    std::ostream* bsnp_diploid_allele_osptr() const { return _bsnp_diploid_allele_osptr.get(); }
    std::ostream* nonref_test_osptr() const { return _nonref_test_osptr.get(); }
    std::ostream* nonref_sites_osptr() const { return _nonref_sites_osptr.get(); }

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
                    std::ostream& os);

private:
    std::ostream& _report_os;
    std::auto_ptr<std::ostream> _counts_osptr;
    std::auto_ptr<std::ostream> _bsnp_diploid_osptr;
    std::auto_ptr<std::ostream> _bsnp_diploid_allele_osptr;
    std::auto_ptr<std::ostream> _nonref_test_osptr;
    std::auto_ptr<std::ostream> _nonref_sites_osptr;
};


#endif

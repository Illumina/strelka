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

/// \file

/// \author Chris Saunders
///
#ifndef __STARLING_STREAMS_BASE_HH
#define __STARLING_STREAMS_BASE_HH

#include "blt_common/blt_streams.hh"
#include "blt_util/bam_dumper.hh"
#include "starling_common/sample_info.hh"
#include "starling_common/starling_shared.hh"
#include "starling_common/starling_types.hh"

#include "boost/shared_ptr.hpp"

#include <vector>


struct starling_streams_base : public blt_streams
{

    typedef blt_streams base_t;

    starling_streams_base(const starling_options& opt,
                          const prog_info& pinfo,
                          const sample_info& si);

    ~starling_streams_base();

    std::ostream*
    bindel_diploid_osptr(const unsigned sample_no) const
    {
        return _bindel_diploid_osptr[sample_no].get();
    }

    std::ostream*
    gvcf_osptr(const unsigned sample_no) const
    {
        return _gvcf_osptr[sample_no];
    }

    bam_dumper*
    realign_bam_ptr(const unsigned sample_no) const
    {
        return _realign_bam_ptr[sample_no].get();
    }

    std::ostream*
    candidate_indel_osptr() const
    {
        return _candidate_indel_osptr.get();
    }

    std::ostream*
    variant_window_osptr(const unsigned window_no) const
    {
        return _window_osptr.at(window_no).get();
    }

protected:
    static
    std::ostream*
    initialize_bindel_file(const starling_options& client_opt,
                           const prog_info& pinfo,
                           const std::string& filename,
                           const char* label=NULL);

    static
    std::ostream*
    initialize_gvcf_file(const starling_options& opt,
                         const prog_info& pinfo,
                         const std::string& filename,
                         const bam_header_t* const header,
                         std::auto_ptr<std::ostream>& os_ptr_auto);

    bam_dumper*
    initialize_realign_bam(const bool is_clobber,
                           const prog_info& pinfo,
                           const std::string& filename,
                           const char* label,
                           const bam_header_t* const header);

    static
    std::ostream*
    initialize_candidate_indel_file(const starling_options& client_opt,
                                    const prog_info& pinfo,
                                    const std::string& filename);

    static
    std::ostream*
    initialize_window_file(const starling_options& opt,
                           const prog_info& pinfo,
                           const avg_window_data& awd,
                           const sample_info& si);

    std::auto_ptr<std::ostream> _bindel_diploid_osptr[MAX_SAMPLE];
    std::ostream* _gvcf_osptr[MAX_SAMPLE];
    std::auto_ptr<std::ostream> _gvcf_osptr_auto[MAX_SAMPLE];
    std::auto_ptr<bam_dumper> _realign_bam_ptr[MAX_SAMPLE];
private:
    std::auto_ptr<std::ostream> _candidate_indel_osptr;
    unsigned _n_samples;

    std::vector<boost::shared_ptr<std::ostream> > _window_osptr;
};


#endif

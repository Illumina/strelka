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

#include "blt_util/bam_dumper.hh"
#include "blt_util/vcf_util.hh"
#include "starling_common/starling_streams_base.hh"
#include "starling_common/gvcf_locus_info.hh"

#include <cassert>

#include <fstream>
#include <iostream>
#include <sstream>



std::ostream*
starling_streams_base::
initialize_bindel_file(const starling_options& opt,
                       const prog_info& pinfo,
                       const std::string& filename,
                       const char* label) {

    const char* const cmdline(opt.cmdline.c_str());

    std::ofstream* fosptr(new std::ofstream);
    std::ofstream& fos(*fosptr);
    open_ofstream(pinfo,filename,"bindel-diploid",opt.is_clobber,fos);

    fos << "# ** " << pinfo.name();
    if(label) fos << " " << label;
    fos << " bindel-diploid file **\n";
    write_file_audit(opt,pinfo,cmdline,fos);
    fos << "#$ INDEL_THETA " << opt.bindel_diploid_theta << "\n";
    fos << "#\n";
    fos << "#$ COLUMNS seq_name pos type ref_upstream ref/indel ref_downstream Q(indel) max_gtype Q(max_gtype) depth alt_reads indel_reads other_reads repeat_unit ref_repeat_count indel_repeat_count\n";

    return fosptr;
}



static
void
write_vcf_filter(std::ostream& os,
                 const char* id,
                 const char* desc) {
    os << "##FILTER=<ID=" << id << ",Description=\"" << desc << "\">\n";
}



static
void
add_gvcf_filters(const gvcf_options& opt,
                 std::ostream& os) {

    using namespace VCF_FILTERS;

    write_vcf_filter(os,get_label(IndelConflict),"Locus is in region with conflicting indel calls.");
    write_vcf_filter(os,get_label(SiteConflict),"Site genotype conflicts with proximal indel call. This is typically a heterozygous SNV call made inside of a heterozygous deletion.");
    
    if(opt.is_max_snv_sb) {
        std::ostringstream oss;
        oss << "SNV strand bias value (SNVSB) exceeds " << opt.max_snv_sb;
        write_vcf_filter(os,get_label(HighSNVSB),oss.str().c_str());
    }
    if(opt.is_max_snv_hpol) {
        std::ostringstream oss;
        oss << "SNV contextual homopolymer length (SNVHPOL) exceeds " <<opt.max_snv_hpol;
        write_vcf_filter(os,get_label(HighSNVHPOL),oss.str().c_str());
    }

    if(opt.is_max_ref_rep) {
        std::ostringstream oss;
        oss << "Indel contains an ellele which occurs in a homopolymer or dinucleotide track with a reference repeat greater than " << opt.max_ref_rep;
        write_vcf_filter(os,get_label(HighRefRep),oss.str().c_str());
    }
}



// return stream, is_malloced_pointer
std::ostream*
starling_streams_base::
initialize_gvcf_file(const starling_options& opt,
                     const prog_info& pinfo,
                     const std::string& filename,
                     std::auto_ptr<std::ostream>& os_ptr_auto) {

    std::ostream* osptr(&std::cout);
    if(filename != "-") {
        std::ofstream* fos_ptr(new std::ofstream);
        open_ofstream(pinfo,filename,"gvcf",opt.is_clobber,*fos_ptr);
        os_ptr_auto.reset(fos_ptr);
        osptr=os_ptr_auto.get();
    }
    std::ostream& os(*osptr);

    const char* const cmdline(opt.cmdline.c_str());

    write_vcf_audit(opt,pinfo,cmdline,os);
    os << "##content=starling small-variant calls\n"
       << "##SnpTheta=" << opt.bsnp_diploid_theta << "\n"
       << "##IndelTheta=" << opt.bindel_diploid_theta << "\n";

    //INFO:
    os << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the region described in this record\">\n";
    os << "##INFO=<ID=BLOCKAVG_min30p3a,Number=0,Type=Flag,Description=\"Non-variant site block. All sites in a block are constrained to be non-variant, have the same filter value, and have all sample values in range [x,y], y <= max(x+3,(x*1.3)). All printed site block sample values are the minimum observed in the region spanned by the block\">\n";

    // site specific:
    os << "##INFO=<ID=SNVSB,Number=1,Type=Float,Description=\"SNV site strand bias\">\n";
    os << "##INFO=<ID=SNVHPOL,Number=1,Type=Integer,Description=\"SNV contextual homopolymer length\">\n";

    // indel specific:
    os << "##INFO=<ID=CIGAR,Number=A,Type=String,Description=\"CIGAR alignment for each alternate indel allele\">\n";
    os << "##INFO=<ID=RU,Number=A,Type=String,Description=\"Smallest repeating sequence unit extended or contracted in the indel allele relative to the reference. RUs are not reported if longer than 20 bases.\">\n";
    os << "##INFO=<ID=REFREP,Number=A,Type=String,Description=\"Number of times RU is repeated in reference.\">\n";
    os << "##INFO=<ID=IREP,Number=A,Type=String,Description=\"Number of times RU is repeated in indel allele.\">\n";

    //FORMAT:
    os << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    os << "##FORMAT=<ID=GQX,Number=1,Type=Integer,Description=\"Minimum of {Genotype quality assuming variant position,Genotype quality assuming non-variant position}\">\n";

    // FILTER:
    add_gvcf_filters(opt.gvcf,os);

    os << vcf_col_label() << "\tFORMAT\tSAMPLE\n";

    return osptr;
}

bam_dumper*
starling_streams_base::
initialize_realign_bam(const bool is_clobber,
                       const prog_info& pinfo,
                       const std::string& filename,
                       const char* label,
                       const bam_header_t* const header) {

    assert(NULL != header);

    // \TODO consider putting extra info into BAM header:
    //
    //fp->header = bam_header_dup((const bam_header_t*)aux);
    //fos << "@PG\tID:" << pinfo.name() << "\tVN:" << pinfo.version() << "\tCL:" << cmdline << "\n";
    
    if(not is_clobber){ // weak clobber test:
        std::ofstream fos;
        open_ofstream(pinfo,filename,label,is_clobber,fos);
    }
    return new bam_dumper(filename.c_str(),header);
}



std::ostream*
starling_streams_base::
initialize_candidate_indel_file(const starling_options& opt,
                                const prog_info& pinfo,
                                const std::string& filename) {

    const char* const cmdline(opt.cmdline.c_str());

    std::ofstream* fosptr(new std::ofstream);
    std::ofstream& fos(*fosptr);
    open_ofstream(pinfo,filename,"candidate-indel",opt.is_clobber,fos);

    fos << "# ** " << pinfo.name();
    fos << " candidate-indel file **\n";
    write_file_audit(opt,pinfo,cmdline,fos);
    fos << "#\n";
    fos << "#$ COLUMNS seq_name pos type length seq\n";

    return fosptr;
}



std::ostream*
starling_streams_base::
initialize_window_file(const starling_options& opt,
                       const prog_info& pinfo,
                       const avg_window_data& awd,
                       const sample_info& si) {

    const char* const cmdline(opt.cmdline.c_str());

    std::ofstream* fosptr(new std::ofstream);
    std::ofstream& fos(*fosptr);
    open_ofstream(pinfo,awd.filename,"variant-window",opt.is_clobber,fos);

    const unsigned fs(awd.flank_size);

    fos << "# ** " << pinfo.name();
    fos << " variant-window file **\n";
    write_file_audit(opt,pinfo,cmdline,fos);
    fos << "#$ FLANK_SIZE " << awd.flank_size << "\n";
    fos << "#\n";
    fos << "#$ COLUMNS seq_name pos";
    static const bool is_tier1(true);
    static const char* win_type[] = {"used","filt","submap"};
    static unsigned n_win_type(sizeof(win_type)/sizeof(char*));
    const unsigned n_samples(si.sample_size());
    for(unsigned s(0);s<n_samples;++s) {
        for(unsigned i(0);i<n_win_type;++i) {
            fos << " " << si.get_prefix(s,is_tier1) << "win" << fs << "_" << win_type[i];
        }
    }
    fos << "\n";

    return fosptr;
}



starling_streams_base::
starling_streams_base(const starling_options& opt,
                      const prog_info& pinfo,
                      const sample_info& si)
    : base_t(opt,pinfo,true)
    , _n_samples(si.sample_size())
    , _window_osptr(opt.variant_windows.size())
{
    assert((_n_samples>0) && (_n_samples<=MAX_SAMPLE));

    for(unsigned i(0);i<_n_samples;++i) _gvcf_osptr[i] = NULL;

    if(opt.is_write_candidate_indels()) {
        _candidate_indel_osptr.reset(initialize_candidate_indel_file(opt,pinfo,opt.candidate_indel_filename));
    }

    const unsigned vs(opt.variant_windows.size());
    for(unsigned i(0);i<vs;++i) {
        _window_osptr[i].reset(initialize_window_file(opt,pinfo,opt.variant_windows[i],si));
    }
}



// dtor here to make auto_ptr work correctly:
starling_streams_base::
~starling_streams_base() {}

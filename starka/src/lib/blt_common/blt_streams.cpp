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

#include "blt_common/blt_streams.hh"
#include "blt_util/digt.hh"
#include "blt_util/vcf_util.hh"

#include <ctime>

#include <fstream>
#include <iostream>
#include <sstream>



static
void
open_ofstream(const prog_info& pinfo,
              const std::string& filename,
              const char* label,
              const bool is_clobber,
              std::ofstream& fos) {

    std::ifstream tis(filename.c_str());
    if(tis && (! is_clobber)) {
        std::ostringstream oss;
        oss << label << " file already exists: " << filename;
        pinfo.usage(oss.str().c_str());
    }

    fos.open(filename.c_str());
    if(!fos) {
        std::ostringstream oss;
        oss << label << " file can't be opened: " << filename;
        pinfo.usage(oss.str().c_str());
    }
}



static
void
write_audit(const blt_options& opt,
            const prog_info& pinfo,
            const char* const cmdline,
            std::ostream& os,
            const char* const prefix = 0){


    if(opt.is_write_variable_metadata){
        if(prefix) os << prefix;
        os << "CMDLINE " << cmdline << "\n";
    }
    if(prefix) os << prefix;
    os << "PROGRAM_VERSION " << pinfo.version() << "\n";
    if(opt.is_write_variable_metadata){
        if(prefix) os << prefix;
        const time_t result(time(0));
        os << "START_TIME " << asctime(localtime(&result));
    }
}



static
void
write_vcf_audit(const blt_options& opt,
                const prog_info& pinfo,
                const char* const cmdline,
                const bam_header_t* const header,
                std::ostream& os) {
    
    const time_t t(time(NULL));

    os << "##fileformat=VCFv4.1\n";
    os << "##fileDate=" << vcf_fileDate << "\n";
    os << "##source=" << pinfo.name() << "\n";
    os << "##source_version=" << pinfo.version() << "\n";
    os << "##startTime=" << asctime(localtime(&t));
    os << "##cmdline=" << cmdline << "\n";
    if(opt.is_samtools_ref_set) {
        os << "##reference=file://" << opt.samtools_ref_seq_file << "\n";
    } else {
        assert(0);
    }

    assert(NULL != header);

    for(int32_t i(0);i<header->n_targets;++i) {
        os << "##contig=<ID=" << header->target_name[i]
           << ",length=" << header->target_len[i] << ">\n";
    }
}



static
void
write_file_audit(const blt_options& opt,
                 const prog_info& pinfo,
                 const char* const cmdline,
                 std::ostream& os) {

    write_audit(opt,pinfo,cmdline,os,"#$ ");
    os << "#\n";
}



static
void
setup_bsnp_file(const blt_options& opt,
                const prog_info& pinfo,
                const char* const cmdline,
                std::auto_ptr<std::ostream>& apo,
                const std::string& filename,
                const char* label,
                const bool is_include_seq_name,
                const char* append_label=NULL) {

    std::ofstream* fosptr(new std::ofstream);
    apo.reset(fosptr);
    std::ofstream& fos(*fosptr);
    open_ofstream(pinfo,filename,label,opt.is_clobber,fos);

    fos << "# ** " << pinfo.name() << " " << label << " file **\n";
    write_file_audit(opt,pinfo,cmdline,fos);
    fos << "#$ SNP_THETA " << opt.bsnp_diploid_theta << "\n";
    if(opt.is_bsnp_diploid_het_bias) {
        fos << "#$ BSNP_DIPLOID_HET_BIAS " << opt.bsnp_diploid_het_bias << "\n";
    }
    fos << "#\n";
    fos << "#$ COLUMNS ";
    if(is_include_seq_name) { fos << "seq_name "; }
    fos << "pos bcalls_used bcalls_filt ref Q(snp) max_gt Q(max_gt) max_gt|poly_site Q(max_gt|poly_site)";

    if(opt.is_print_used_allele_counts) {
        for(unsigned b(0);b<N_BASE;++b){
            fos << " " << id_to_base(b) << "_used";
        }
    }

    fos << " snv_sb snv_hpol";
    
    if(opt.is_print_all_poly_gt) {
        for(unsigned gt(0);gt<DIGT::SIZE;++gt){
            fos << " P(" << DIGT::label(gt) << "|ps)";
        }
    }

    if(NULL != append_label) {
        fos << " " << append_label;
    }

    fos << "\n";
}



static
void
setup_nonref_output(const blt_options& opt,
                    const prog_info& pinfo,
                    std::auto_ptr<std::ostream>& osptr,
                    const char* filename,
                    const char* label) {

    const char* const cmdline(opt.cmdline.c_str());

    std::ofstream* fosptr(new std::ofstream);
    osptr.reset(fosptr);
    std::ofstream& fos(*fosptr);
    open_ofstream(pinfo,filename,label,opt.is_clobber,fos);
    
    fos << "# ** " << pinfo.name() << " nonref allele test file **\n";
    write_file_audit(opt,pinfo,cmdline,fos);
    fos << "#$ COLUMNS seq_name pos bcalls_used bcalls_filt ref Q(snv) max_gt Q(max_gt)";
    
    //        if(opt.is_print_used_allele_counts) {
    for(unsigned b(0);b<N_BASE;++b){
        fos << " " << id_to_base(b) << "_used";
    }
    
    for(unsigned b(0);b<N_BASE;++b){
        fos << " " << id_to_base(b) << "_meanQ";
    }
    //}
    
    fos << "\n";
}



void
blt_streams::
write_file_audit(const blt_options& opt,
                 const prog_info& pinfo,
                 const char* const cmdline,
                 std::ostream& os) {

    ::write_file_audit(opt,pinfo,cmdline,os);
}



void
blt_streams::
write_vcf_audit(const blt_options& opt,
                const prog_info& pinfo,
                const char* const cmdline,
                const bam_header_t* const header,
                std::ostream& os) {

    ::write_vcf_audit(opt,pinfo,cmdline,header,os);
}



void
blt_streams::
open_ofstream(const prog_info& pinfo,
              const std::string& filename,
              const char* label,
              const bool is_clobber,
              std::ofstream& fos) {

    ::open_ofstream(pinfo,filename,label,is_clobber,fos);
}



blt_streams::
blt_streams(const blt_options& opt,
            const prog_info& pinfo,
            const bool is_include_seq_name) {

    const char* const cmdline(opt.cmdline.c_str());

    if(! opt.report_filename.empty()) {
        std::ofstream* fosptr(new std::ofstream);
        _report_osptr.reset(fosptr);
        std::ofstream& fos(*fosptr);
        open_ofstream(pinfo,opt.report_filename,"report",opt.is_clobber,fos);

        write_audit(opt,pinfo,cmdline,fos);
    }

    if(opt.is_counts) {
        std::ofstream* fosptr(new std::ofstream);
        _counts_osptr.reset(fosptr);
        std::ofstream& fos(*fosptr);
        open_ofstream(pinfo,opt.counts_filename,"counts",opt.is_clobber,fos);


        fos << "# ** " << pinfo.name() << " counts file **\n";
        write_file_audit(opt,pinfo,cmdline,fos);
        fos << "#$ COLUMNS pos A_used C_used G_used T_used unused\n";
    }

    if(opt.is_nonref_test()){
        setup_nonref_output(opt,pinfo,_nonref_test_osptr,opt.nonref_test_filename.c_str(),"nonref test");
    }

    if(opt.is_nonref_sites()){
        setup_nonref_output(opt,pinfo,_nonref_sites_osptr,opt.nonref_sites_filename.c_str(),"nonref sites");
    }

    if(opt.is_nonref_sites()){
        std::ofstream* fosptr(new std::ofstream);
        _nonref_test_osptr.reset(fosptr);
        std::ofstream& fos(*fosptr);
        open_ofstream(pinfo,opt.nonref_test_filename,"nonref test",opt.is_clobber,fos);

        fos << "# ** " << pinfo.name() << " nonref allele test file **\n";
        write_file_audit(opt,pinfo,cmdline,fos);
        fos << "#$ COLUMNS seq_name pos bcalls_used bcalls_filt ref Q(snv) max_gt Q(max_gt)";

        //        if(opt.is_print_used_allele_counts) {
        for(unsigned b(0);b<N_BASE;++b){
            fos << " " << id_to_base(b) << "_used";
        }

        for(unsigned b(0);b<N_BASE;++b){
            fos << " " << id_to_base(b) << "_meanQ";
        }
        //}

        fos << "\n";
    }

    if(opt.is_bsnp_diploid_file){
        setup_bsnp_file(opt,pinfo,cmdline,_bsnp_diploid_osptr,
                        opt.bsnp_diploid_filename,"bsnp-diploid",is_include_seq_name);          
    }

    if(opt.is_bsnp_diploid_allele_file){
        setup_bsnp_file(opt,pinfo,cmdline,_bsnp_diploid_allele_osptr,
                        opt.bsnp_diploid_allele_filename,"bsnp-diploid sites",is_include_seq_name);
    }
}



// this should remain in the cpp so that auto_ptr works correctly
blt_streams::
~blt_streams() { }

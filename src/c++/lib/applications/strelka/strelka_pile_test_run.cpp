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
/// note coding convention for all ranges '_pos fields' is:
/// XXX_begin_pos is zero-indexed position at the begining of the range
/// XXX_end_pos is zero-index position 1 step after the end of the range
///


#include "position_somatic_snv_strand_grid.hh"
#include "strelka_sample_type.hh"
#include "strelka_pile_test_run.hh"

#include "blt_util/istream_line_splitter.hh"
#include "blt_util/parse_util.hh"
#include "blt_util/seq_util.hh"

#include <cassert>
#include <cctype>

#include <iostream>



strelka_pile_caller::
strelka_pile_caller(strelka_options& opt,
                    std::ostream& os)
    : _opt(opt)
    , _os(os)
{
    // set default parameters:
    opt.is_bsnp_diploid_file = true;
    opt.bsnp_diploid_theta = 0.001;
    opt.somatic_snv_rate       = 0.000001;
    opt.shared_site_error_rate = 0.0000005;
    opt.shared_site_error_strand_bias_fraction = 0.5;

    static const std::string ref_seq("ACGT");
    reference_contig_segment ref;
    ref.seq()=ref_seq;
    _dopt_ptr.reset(new strelka_deriv_options(opt,ref));
}



void
strelka_pile_caller::
call(
    const bool is_somatic_gvcf,
    const unsigned pos,
    snp_pos_info& norm_pi,
    snp_pos_info& tumor_pi) {

    static dependent_prob_cache dpcache;

    const char ref_base(norm_pi.ref_base);

    // recreate data caches:
    extra_position_data norm_epd;
    extra_position_data tumor_epd;

    static const bool is_dep(false);
    const bool is_include_tier2(false);
    extended_pos_data normald(&norm_pi,norm_epd,
                              ref_base,_opt,dpcache,is_dep,is_include_tier2);
    extended_pos_data tumord(&tumor_pi,tumor_epd,
                             ref_base,_opt,dpcache,is_dep,is_include_tier2);

    //    somatic_snv_genotype sgt;
    somatic_snv_genotype_grid sgtg;
    _dopt_ptr->sscaller_strand_grid().position_somatic_snv_call(normald.good_epi,
                                                                tumord.good_epi,
                                                                NULL,
                                                                NULL,
                                                                is_somatic_gvcf,
                                                                sgtg);

    if (! (sgtg.is_output() || is_somatic_gvcf)) return;

    static const char chrom_name[] = "sim";
    _os << chrom_name << '\t'
        << pos << '\t'
        << ".";

    write_vcf_somatic_snv_genotype_strand_grid(_opt,sgtg,is_somatic_gvcf,
                                               normald,
                                               tumord,
                                               normald,
                                               tumord,
                                               _os);

    _os << "\n";
}



static
void
load_pi(const char ref_base,
        const char* read,
        const uint8_t* qual,
        snp_pos_info& pi) {

    pi.clear();
    pi.ref_base=ref_base;

    const unsigned len(strlen(read));
    for (unsigned i(0); i<len; ++i) {
        const bool is_fwd(isupper(read[i]));
        const uint8_t base_id(base_to_id(toupper(read[i])));
        assert(qual[i]>=33);
        pi.calls.push_back(base_call(base_id,qual[i]-33,is_fwd,
                                     1,1,false,false));
    }
}


// call sites from an external simulator:
void
strelka_pile_test_run(strelka_options& opt) {

    strelka_pile_caller pcall(opt,std::cout);

    istream_line_splitter dparse(std::cin);

    snp_pos_info norm_pi;
    snp_pos_info tumor_pi;

    while (dparse.parse_line()) {

        assert(6 == dparse.n_word());

        const char* pcopy(dparse.word[0]);
        const unsigned pos(illumina::blt_util::parse_unsigned(pcopy));
        const char ref_base(*dparse.word[1]);

        {
            const char* normbase(dparse.word[2]);
            const uint8_t* normqual((uint8_t*) dparse.word[3]);
            load_pi(ref_base,normbase,normqual,norm_pi);
        }

        {
            const char* tumorbase(dparse.word[4]);
            const uint8_t* tumorqual((uint8_t*) dparse.word[5]);
            load_pi(ref_base,tumorbase,tumorqual,tumor_pi);
        }

        static const bool is_somatic_gvcf(false);
        pcall.call(is_somatic_gvcf, pos,norm_pi,tumor_pi);
    }
}

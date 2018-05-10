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

///
/// \author Chris Saunders
///
/// note coding convention for all ranges '_pos fields' is:
/// XXX_begin_pos is zero-indexed position at the beginning of the range
/// XXX_end_pos is zero-index position 1 step after the end of the range
///


#include "position_somatic_snv_strand_grid.hh"
#include "position_somatic_snv_strand_grid_vcf.hh"
#include "strelka_pile_test_run.hh"

#include "blt_util/istream_line_splitter.hh"
#include "blt_util/parse_util.hh"
#include "blt_util/seq_util.hh"
#include "starling_common/PileupCleaner.hh"

#include <cassert>
#include <cctype>

#include <iostream>

#include "strelka_common/StrelkaSampleSetSummary.hh"



strelka_pile_caller::
strelka_pile_caller(strelka_options& opt,
                    std::ostream& os)
    : _opt(opt)
    , _os(os)
{
    // set default parameters:
    opt.bsnp_diploid_theta = 0.001;
    opt.somatic_snv_rate       = 0.000001;
    opt.shared_site_error_rate = 0.0000005;
//    opt.shared_site_error_rate = 0.0000000005;
//    opt.shared_site_error_strand_bias_fraction = 0.5;
    opt.shared_site_error_strand_bias_fraction = 0.0;

    static const std::string ref_seq("ACGT");
    reference_contig_segment ref;
    ref.seq()=ref_seq;
    _dopt_ptr.reset(new strelka_deriv_options(opt));
}



void
strelka_pile_caller::
call(
    const bool is_somatic_gvcf,
    const unsigned pos,
    snp_pos_info& norm_pi,
    snp_pos_info& tumor_pi)
{
    static PileupCleaner pileupCleaner(_opt);

    // recreate data caches:
    CleanedPileup norm_cpi;
    CleanedPileup tumor_cpi;

    const bool is_include_tier2(false);
    pileupCleaner.CleanPileup(norm_pi,is_include_tier2,norm_cpi);
    pileupCleaner.CleanPileup(tumor_pi,is_include_tier2,tumor_cpi);

    //    somatic_snv_genotype sgt;
    somatic_snv_genotype_grid sgtg;
    _dopt_ptr->sscaller_strand_grid().position_somatic_snv_call(norm_cpi.getExtendedPosInfo(),
                                                                tumor_cpi.getExtendedPosInfo(),
                                                                nullptr,
                                                                nullptr,
                                                                is_somatic_gvcf,
                                                                sgtg);

    if (! (sgtg.is_output() || is_somatic_gvcf)) return;

    static const char chrom_name[] = "sim";
    _os << chrom_name << '\t'
        << pos << '\t'
        << ".";

    write_vcf_somatic_snv_genotype_strand_grid(_opt, *(_dopt_ptr), sgtg, is_somatic_gvcf, norm_cpi,
                                               tumor_cpi, norm_cpi, tumor_cpi, 0, 0, _os);

    _os << "\n";
}



static
void
load_pi(const char ref_base,
        const char* read,
        const uint8_t* qual,
        snp_pos_info& pi)
{
    pi.clear();
    pi.set_ref_base(ref_base);

    const unsigned len(strlen(read));
    for (unsigned i(0); i<len; ++i)
    {
        const bool is_fwd(isupper(read[i])!=0);
        const uint8_t base_id(base_to_id(toupper(read[i])));
        assert(qual[i]>=33);
        pi.calls.push_back(base_call(base_id,qual[i]-33,is_fwd,
                                     1,1,false,false,false));
    }
}


// call sites from an external simulator:
void
strelka_pile_test_run(
    strelka_options& opt)
{
    strelka_pile_caller pcall(opt,std::cout);

    istream_line_splitter dparse(std::cin);

    snp_pos_info norm_pi;
    snp_pos_info tumor_pi;

    while (dparse.parse_line())
    {

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

// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

/// \file
///
/// \author Chris Saunders
///
/// note coding convention for all ranges '_pos fields' is:
/// XXX_begin_pos is zero-indexed position at the begining of the range
/// XXX_end_pos is zero-index position 1 step after the end of the range
///


#include "position_somatic_snv_strand_grid.hh"
#include "strelka_sample_type.hh"

#include "blt_util/istream_line_splitter.hh"
#include "blt_util/qscore.hh"
#include "blt_util/seq_util.hh"
#include "strelka/strelka_pile_test_run.hh"
#include "strelka/strelka_sim_test.hh"

#include "boost/lexical_cast.hpp"
#include "boost/random.hpp"
#include "boost/random/binomial_distribution.hpp"
#include "boost/random/variate_generator.hpp"

#include <cassert>
#include <cctype>

#include <fstream>
#include <iostream>


typedef boost::mt19937 gen_t;
typedef boost::poisson_distribution<> dist_t;
typedef boost::variate_generator<gen_t&,dist_t> vgen_t;


// shared rng:
gen_t gen;//(static_cast<unsigned> (std::time(0));

// shared [0,1]
typedef boost::uniform_real<> rdist_t;
rdist_t rdist;
boost::variate_generator<gen_t&,rdist_t> uran(gen,rdist);

// shared 3range:
typedef boost::uniform_smallint<> sudist_t;
sudist_t dist3i(0,2);
boost::variate_generator<gen_t&,sudist_t> ran3i(gen,dist3i);



static
unsigned
get_binom(const unsigned n,
          const float p) {
    typedef boost::binomial_distribution<int,float> bdist_t;
    typedef boost::variate_generator<gen_t&,bdist_t> bvgen_t;

    bdist_t bdist(n,p);
    bvgen_t bvgen(gen,bdist);
    return static_cast<unsigned>(bvgen());
}



static
inline
uint8_t
get_mut_base_id(const uint8_t base_id) {
    uint8_t id(ran3i());
    if (id>=base_id) id += 1;
    return id;
}



static
uint8_t
get_obs_base_id(const uint8_t true_id,
                const uint8_t qval) {

    if (uran() >= qphred_to_error_prob(qval)) {
        return true_id;
    }

    return get_mut_base_id(true_id);
}



/// \brief pull a random variate from a discrete cdf
///
static
unsigned
random_cdf_variate(const double cdf[],
                   const unsigned N) {

    const double* lbp(std::lower_bound(cdf,cdf+N,uran()));
    return std::min(static_cast<unsigned>(lbp-cdf),N-1);
}



struct qval_distro {

    qval_distro(
        const uint8_t constval = 30)
        : _is_const(true),
          _constval(constval),
          _qsize(0)
    {}

    qval_distro(const char* distro_file)
        : _is_const(false),
          _constval(0),
          _qsize(0)
    {
        double total(0);

        assert(distro_file);
        std::ifstream ifs(distro_file);
        if (! ifs) {
            std::cerr << "ERROR: can't open file '" << distro_file << "'\n";
            exit(EXIT_FAILURE);
        }

        istream_line_splitter dparse(ifs);

        while (dparse.parse_line()) {
            assert(dparse.n_word() > 0);
            const char* pcopy(dparse.word[0]);
            if (strlen(pcopy) && pcopy[0] == '#') continue;

            assert(_qsize+1 < MAX_QVAL);
            assert(dparse.n_word() > 1);

            char* word0(dparse.word[0]);
            char* word1(dparse.word[1]);
            unsigned id(boost::lexical_cast<unsigned>(word0));
            assert(id<MAX_QVAL);
            _qval_id[_qsize] = id;
            _qval_cdf[_qsize] = boost::lexical_cast<unsigned>(word1);
            total += _qval_cdf[_qsize];
            _qsize++;
        }

        for (unsigned i(0); i<_qsize; ++i) {
            _qval_cdf[i] /= total;
            if (i) _qval_cdf[i] += _qval_cdf[i-1];
        }
    }

    uint8_t
    get() const
    {
        if (_is_const)
        {
            return _constval;
        }
        else
        {
            return _qval_id[random_cdf_variate(_qval_cdf,_qsize)];
        }
    }

private:
    bool _is_const;
    uint8_t _constval;

    enum { MAX_QVAL=80 };
    unsigned _qsize;
    uint8_t _qval_id[MAX_QVAL];
    double _qval_cdf[MAX_QVAL];
};



static
void
sim_sample_pi(vgen_t& cov_gen,
              const qval_distro& qdist,
              const uint8_t ref_id,
              const uint8_t alt_id,
              const float alt_freq,
              snp_pos_info& pi) {

    const unsigned all_cov(cov_gen());
    const unsigned fwd_cov(get_binom(all_cov,0.5));
    const unsigned rev_cov(all_cov-fwd_cov);

    unsigned fwd_alt(0);
    unsigned rev_alt(0);
    if (alt_freq > 0) {
        fwd_alt=get_binom(fwd_cov,alt_freq);
        rev_alt=get_binom(rev_cov,alt_freq);
    }

    pi.clear();
    pi.ref_base=id_to_base(ref_id);

    for (unsigned i(0); i<all_cov; ++i) {
        const uint8_t qval(qdist.get());

        uint8_t true_id(ref_id);
        if ( (i < fwd_alt) || ((i >= fwd_cov) && (i < (fwd_cov+rev_alt)))) {
            true_id=alt_id;
        }

        const uint8_t obs_id(get_obs_base_id(true_id,qval));
        const bool is_fwd(i < fwd_cov);
        pi.calls.push_back(base_call(obs_id,qval,is_fwd,
                                     1,1,false,false));
    }
}



// completely self-contained version of pile simulator:
//
void
strelka_site_sim(strelka_options& opt,
                 strelka_site_sim_options& sim_opt) {

    gen.seed(sim_opt.seed);

    snp_pos_info norm_pi;
    snp_pos_info tumor_pi;

    std::auto_ptr<qval_distro> qdistptr;
    if (sim_opt.qval_file.empty())
    {
        qdistptr.reset(new qval_distro());
    }
    else
    {
        qdistptr.reset(new qval_distro(sim_opt.qval_file.c_str()));
    }

    dist_t ndist(sim_opt.ncov);
    dist_t tdist(sim_opt.tcov);
    vgen_t ngen(gen,ndist);
    vgen_t tgen(gen,tdist);

    static const uint8_t ref_id(0);

    bool is_ofs(false);
    std::ofstream ofs;
    if (! sim_opt.oracle_file.empty()) {
        is_ofs=true;
        ofs.open(sim_opt.oracle_file.c_str());
    }

    strelka_pile_caller scall(opt,std::cout);

    for (unsigned i(0); i<sim_opt.total_sites; ++i) {

        const unsigned pos(i+1);

        uint8_t nalt_id(ref_id);
        uint8_t talt_id(ref_id);
        float nalt_freq(0.);
        float talt_freq(0.);

        // test for alternate states
        if (sim_opt.mode == SIM_RANDOM) {
            if (uran() < opt.shared_site_error_rate) {
                sim_opt.mode=SIM_NOISE;
            } else if (uran() < sim_opt.ssnv_prior) {
                sim_opt.mode=SIM_SOMATIC;
            } else if (uran() < opt.bsnp_diploid_theta) {
                sim_opt.mode=SIM_GERMLINE;
            } else {
                sim_opt.mode=SIM_REF;
            }
        }

        //bool is_nonref(false);
        if (sim_opt.mode == SIM_NOISE) {
            nalt_id=get_mut_base_id(nalt_id);
            talt_id=nalt_id;
            nalt_freq=uran();
            talt_freq=nalt_freq;

            if (is_ofs) ofs << pos << "\tNOISE\t" << nalt_freq << "\n";
            //is_nonref=true;

        } else if (sim_opt.mode == SIM_SOMATIC) {
            talt_id=get_mut_base_id(nalt_id);
            talt_freq=0.5*sim_opt.tumor_purity;

            if (is_ofs) ofs << pos << "\tSOMATIC\t" << talt_freq << "\n";
            //is_nonref=true;

        } else if (sim_opt.mode == SIM_GERMLINE) {
            nalt_id=get_mut_base_id(nalt_id);
            talt_id=nalt_id;

            static const double one_third(1./3.);
            if (uran() <= one_third) { nalt_freq=1.; }
            else                    { nalt_freq=0.5; }
            talt_freq=nalt_freq;

            if (is_ofs) ofs << pos << "\tGERMLINE\t" << nalt_freq << "\n";
            //is_nonref=true;
        }

        sim_sample_pi(ngen,*qdistptr,ref_id,nalt_id,nalt_freq,norm_pi);
        sim_sample_pi(tgen,*qdistptr,ref_id,talt_id,talt_freq,tumor_pi);

        scall.call(sim_opt.is_somatic_gvcf,i+1,norm_pi,tumor_pi);
    }
}


#if 0
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



void
strelka_pile_test_run(strelka_options& opt) {

    static const char chrom_name[] = "sim";

    // set default parameters:
    opt.is_bsnp_diploid = true;
    opt.bsnp_diploid_theta = 0.001;
    opt.somatic_snv_rate       = 0.000001;
    opt.shared_site_error_rate = 0.0000005;

    const std::string ref_seq("ACGT");

    const strelka_deriv_options dopt(opt,ref_seq);

    istream_line_splitter dparse(std::cin);

    snp_pos_info norm_pi;
    snp_pos_info tumor_pi;

    // recreate data caches:
    extra_position_data norm_epd;
    extra_position_data tumor_epd;
    extra_position_data tier2_epd[MAX_SAMPLE];

    static dependent_prob_cache dpcache;

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

        static const unsigned n_tier(2);

        std::auto_ptr<extended_pos_data> normald_ptr[n_tier];
        std::auto_ptr<extended_pos_data> tumord_ptr[n_tier];

        extra_position_data* normal_epd_ptr[n_tier] = { &(norm_epd) , &(tier2_epd[STRELKA_SAMPLE_TYPE::NORMAL]) };
        extra_position_data* tumor_epd_ptr[n_tier] = { &(tumor_epd) , &(tier2_epd[STRELKA_SAMPLE_TYPE::TUMOR]) };

        for (unsigned t(0); t<n_tier; ++t) {
            static const bool is_dep(false);
            const bool is_include_tier2(t!=0);
            if (is_include_tier2) continue;
            normald_ptr[t].reset(new extended_pos_data(&norm_pi,*(normal_epd_ptr[t]),
                                                       ref_base,opt,dpcache,is_dep,is_include_tier2));
            tumord_ptr[t].reset(new extended_pos_data(&tumor_pi,*(tumor_epd_ptr[t]),
                                                      ref_base,opt,dpcache,is_dep,is_include_tier2));
        }

        //    somatic_snv_genotype sgt;
        somatic_snv_genotype_grid sgtg;

        static const extended_pos_info* normal_epi_t2_ptr(NULL);
        static const extended_pos_info* tumor_epi_t2_ptr(NULL);

        dopt.sscaller_strand_grid().position_somatic_snv_call(normald_ptr[0]->good_epi,
                                                              tumord_ptr[0]->good_epi,
                                                              normal_epi_t2_ptr,
                                                              tumor_epi_t2_ptr,
                                                              sgtg);

        if (! sgtg.is_snv) continue;

        std::ostream& bos(std::cout);

        // OUTPUT_VCF:
        bos << chrom_name << '\t'
            << pos << '\t'
            << ".";

        write_vcf_somatic_snv_genotype_strand_grid(opt,sgtg,
                                                   *(normald_ptr[0]),
                                                   *(tumord_ptr[0]),
                                                   *(normald_ptr[0]),
                                                   *(tumord_ptr[0]),
                                                   bos);

        bos << "\n";
    }
}
#endif

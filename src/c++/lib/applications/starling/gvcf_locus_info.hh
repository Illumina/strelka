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


#include "blt_common/position_snp_call_pprob_digt.hh"
#include "blt_util/align_path.hh"
#include "starling_common/starling_indel_call_pprob_digt.hh"

#include <bitset>
#include <iosfwd>
#include <map>


namespace VCF_FILTERS
{

enum index_t
{
    IndelConflict,
    SiteConflict,
    PloidyConflict,
    OffTarget,
    LowGQX,
    LowQscoreHetSNP,
    LowQscoreHomSNP,
    LowQscoreHetAltSNP,
    LowQscoreHetIns,
    LowQscoreHomIns,
    LowQscoreHetAltIns,
    LowQscoreHetDel,
    LowQscoreHomDel,
    LowQscoreHetAltDel,
    PhasingConflict,
    HighBaseFilt,
    HighDepth,
    HighSNVSB,
    HighSNVHPOL,
    HighRefRep,
    SIZE
};

inline
const char*
get_label(const unsigned idx)
{
    switch (idx)
    {
    case HighDepth:
        return "HighDepth";
    case LowGQX:
        return "LowGQX";
    case LowQscoreHetSNP:
        return "LowGQXHetSNP";
    case LowQscoreHomSNP:
        return "LowGQXHomSNP";
    case LowQscoreHetAltSNP:
        return "LowGQXHetAltSNP";
    case LowQscoreHetIns:
        return "LowGQXHetIns";
    case LowQscoreHomIns:
        return "LowGQXHomIns";
    case LowQscoreHetAltIns:
        return "LowGQXHetAltIns";
    case LowQscoreHetDel:
        return "LowGQXHetDel";
    case LowQscoreHomDel:
        return "LowGQXHomDel";
    case LowQscoreHetAltDel:
        return "LowGQXHetAltDel";
    case PhasingConflict:
        return "PhasingConflict";
    case HighSNVSB:
        return "HighSNVSB";
    case HighSNVHPOL:
        return "HighSNVHPOL";
    case HighBaseFilt:
        return "HighDPFRatio";
    case HighRefRep:
        return "HighREFREP";
    case IndelConflict:
        return "IndelConflict";
    case SiteConflict:
        return "SiteConflict";
    case PloidyConflict:
        return "PLOIDY_CONFLICT";
    case OffTarget:
        return "OffTarget";
    default:
        assert(false && "Unknown VCF filter value");
        return nullptr;
    }
}
}



struct shared_call_info
{
    shared_call_info()
    {
        clear();
    }
    virtual ~shared_call_info() {}

    void
    set_filter(const VCF_FILTERS::index_t i)
    {
        filters.set(i);
    }

    void
    unset_filter(const VCF_FILTERS::index_t i)
    {
        filters.reset(i);
    }

    void
    write_filters(std::ostream& os) const;

    void
    clear()
    {
        filters.reset();
        gq = gqx = 0;
    }

    int gqx=0;
    int gq=0;
    double strand_bias = 0;

    std::bitset<VCF_FILTERS::SIZE> filters;
};


std::ostream& operator<<(std::ostream& os,const shared_call_info& shmod);

struct shared_indel_call_info : public shared_call_info
{
    shared_indel_call_info(const indel_key& ik,
                           const indel_data& id,
                           const starling_indel_report_info iri,
                           const starling_indel_sample_report_info& isri)
        : _ik(ik)
        , _id(id)
        , _iri(iri)
        , _isri(isri)
    {
    }
    const indel_key _ik;
    const indel_data _id;
    // TODO: make the indel overlapping code create a new call, then revert this to const
    starling_indel_report_info _iri;
    const starling_indel_sample_report_info _isri;

    void set_hap_cigar(
        const unsigned lead=1,
        const unsigned trail=0);

    ALIGNPATH::path_t cigar;
};



struct digt_indel_call : public shared_indel_call_info
{
    digt_indel_call(const indel_key& ik,
                    const indel_data& id,
                    const starling_indel_report_info& iri,
                    const starling_indel_sample_report_info& isri,
                    const starling_diploid_indel_core& dindel)
        : shared_indel_call_info(ik, id, iri, isri)
        , _dindel(dindel)
    {
    }

    // TODO: Make indel_overlapper create new call objects, then revert this to const
    starling_diploid_indel_core _dindel;

    // The empirically calibrated quality-score of an indel, if -1 no q-score has been reported
    int Qscore = -1;

    unsigned max_gt=0;
};



namespace MODIFIED_SITE_GT
{

enum index_t
{
    NONE,
    UNKNOWN,
    ZERO,
    ONE
};

inline
const char*
get_label(const unsigned idx)
{
    switch (static_cast<index_t>(idx))
    {
    case ZERO:
        return "0";
    case ONE:
        return "1";
    case UNKNOWN:
        return ".";
    default:
        assert(false && "Unknown site GT value");
        return nullptr;
    }
}
}

struct digt_call_info : public shared_call_info
{
    digt_call_info()
    {
        clear();
    }

    void
    clear()
    {
        shared_call_info::clear();
        is_unknown=true;
        is_covered=false;
        is_used_covered=false;
        is_zero_ploidy=false;
        is_phased_region=false;
        is_phasing_insufficient_depth=false;
        modified_gt=MODIFIED_SITE_GT::NONE;
        Qscore = -1;
    }

    bool
    is_gqx() const
    {
        return ((!is_unknown) && is_used_covered && (!is_zero_ploidy));
    }

    bool is_unknown;
    bool is_covered;
    bool is_used_covered;
    bool is_zero_ploidy;
    bool is_phased_region;
    bool is_phasing_insufficient_depth;

    MODIFIED_SITE_GT::index_t modified_gt;
    unsigned max_gt;


    // The empirically calibrated quality-score of the site, if -1 not q-score has been reported
    int Qscore = -1;
};

struct continuous_site_call : public shared_call_info
{
    continuous_site_call(int totalDepth, int alleleDepth, BASE_ID::index_t base)
        : _totalDepth(totalDepth)
        , _alleleDepth(alleleDepth)
        , _base(base)
    {
    }

    double variant_frequency() const
    {
        return _totalDepth > 0 ? _alleleDepth / (double)_totalDepth : 0.0;
    }
    int _totalDepth;
    int _alleleDepth;
    BASE_ID::index_t _base;
};

struct continuous_indel_call : public shared_indel_call_info
{
    continuous_indel_call(unsigned totalDepth, unsigned alleleDepth,
                          const indel_key& ik, const indel_data& id, const starling_indel_report_info& iri, const starling_indel_sample_report_info& isri)
        : shared_indel_call_info(ik, id, iri, isri)
        , _totalDepth(totalDepth)
        , _alleleDepth(alleleDepth)
    {
        set_hap_cigar(0,0);
    }

    double variant_frequency() const
    {
        return _totalDepth > 0 ? _alleleDepth / (double)_totalDepth : 0.0;
    }
    unsigned _totalDepth;
    unsigned _alleleDepth;
};



std::ostream& operator<<(std::ostream& os,const digt_call_info& smod);

struct indel_info
{
    explicit indel_info(const pos_t init_pos)
        : pos(init_pos)
    {
    }

    virtual bool is_forced_output() const = 0;
    virtual bool is_indel() const = 0;
    virtual void set_filter(VCF_FILTERS::index_t filter) = 0;
    pos_t pos;
    // the EXCLUSIVE end of the variant (i.e. open)
    virtual pos_t end() const = 0;
};


struct digt_indel_info : public indel_info
{
    digt_indel_info(const pos_t init_pos,
                    const indel_key& init_ik,
                    const indel_data& init_id,
                    const starling_diploid_indel_core& init_dindel,
                    const starling_indel_report_info& init_iri,
                    const starling_indel_sample_report_info& init_isri) : indel_info(init_pos)
    {
        _calls.emplace_back(init_ik, init_id, init_iri, init_isri, init_dindel);
    }

    bool is_forced_output() const override
    {
        return std::any_of(_calls.begin(), _calls.end(),
                           [](const digt_indel_call& x)
        {
            return x._dindel.is_forced_output;
        });
    }
    bool is_indel() const override
    {
        return std::any_of(_calls.begin(), _calls.end(), [](const digt_indel_call& x)
        {
            return x._dindel.is_indel;
        });
    }
    void set_filter(VCF_FILTERS::index_t filter) override
    {
        for (auto& x : _calls) x.set_filter(filter);
    }
    pos_t end() const override;

    void add_overlap(const reference_contig_segment& ref, digt_indel_info& overlap);

    const char*
    get_gt() const
    {
        if (this->is_hetalt())
        {
            return "1/2";
        }
        else if (first()._dindel.is_haploid())
        {
            using namespace STAR_DIINDEL;

            switch (first().max_gt)
            {
            case NOINDEL:
                return "0";
            case HOM:
                return "1";
            default:
                assert(false && "Invalid indel genotype index");
                return "X";
            }
        }
        return STAR_DIINDEL::get_gt_label(first().max_gt);
    }

    void
    dump_pl() const
    {

    }

    bool
    is_hetalt() const
    {
        return (_is_overlap);
    }

    bool
    is_het() const
    {
        return (static_cast<int>(_calls.front().max_gt)>1);
    }

    // the site ploidy within the indel at offset x
    unsigned
    get_ploidy(const unsigned offset) const
    {
        auto& dindel(first()._dindel);
        if (dindel.is_noploid()) return 0;

        if (!is_hetalt())
        {
            using namespace STAR_DIINDEL;
            switch (dindel.max_gt)
            {
            case HOM:
                return 0;
            case HET:
                return 1;
            case NOINDEL:
                return (dindel.is_haploid() ? 1 : 2);
            }
            assert(0);
        }
        else
        {
            assert(offset<ploidy.size());
            return ploidy[offset];
        }
        return 2;
    }
    std::map<std::string, double>
    get_indel_qscore_features(const double chrom_depth) const;

    /// represent site ploidy over the reference span of the overlapping indel set in the event of overlap:
    std::vector<unsigned> ploidy;

    // used to flag hetalt
    bool _is_overlap=false;

    std::vector<digt_indel_call> _calls;

    const digt_indel_call& first() const
    {
        return _calls.front();
    }
    digt_indel_call& first()
    {
        return _calls.front();
    }
};



//Data structure defining parameters for a single site to be used for writing in gvcf_aggregator
struct site_info
{
    site_info(const pos_t init_pos,
              const char init_ref,
              const snp_pos_info& good_pi,
              const int used_allele_count_min_qscore,
              const bool is_forced_output = false)
    {
        pos=(init_pos);
        ref=(init_ref);
        forcedOutput = is_forced_output;
        good_pi.get_known_counts(known_counts,used_allele_count_min_qscore);
        spanning_deletions = good_pi.n_spandel;
    }

    site_info() {}

    virtual bool is_snp() const = 0;
    virtual void set_filter(VCF_FILTERS::index_t filter) = 0;
    virtual bool is_nonref() const = 0;


    pos_t pos = 0;
    char ref = 'N';
    std::array<unsigned,N_BASE> known_counts = {{}};
    unsigned n_used_calls = 0;
    unsigned n_unused_calls = 0;
    unsigned hpol = 0;


    unsigned spanning_deletions;
    bool Unphasable = false;        // Set to true if the site should never be included in a phasing block
    bool forcedOutput = false;
};

struct digt_site_info : public site_info
{
    digt_site_info(const pos_t init_pos,
                   const char init_ref,
                   const snp_pos_info& good_pi,
                   const int used_allele_count_min_qscore,
                   const bool is_forced_output = false)
        : site_info(init_pos, init_ref, good_pi, used_allele_count_min_qscore, is_forced_output)
    {}

    digt_site_info() {}

    bool is_snp() const override
    {
        return dgt.is_snp;
    }
    void set_filter (VCF_FILTERS::index_t filter) override
    {
        smod.set_filter(filter);
    }


    const char*
    get_gt() const
    {
        if       (smod.modified_gt != MODIFIED_SITE_GT::NONE)
        {
            return MODIFIED_SITE_GT::get_label(smod.modified_gt);
        }
        else if (is_print_unknowngt())
        {
            return ".";
        }
        else
        {
            const unsigned print_gt(smod.max_gt);
            return DIGT::get_vcf_gt(print_gt,dgt.ref_gt);
        }
    }

    std::map<std::string, double>
    get_site_qscore_features(const double chrom_depth) const;

    bool
    is_het() const
    {
        unsigned print_gt(smod.max_gt);
        return DIGT::is_het(print_gt);
    }

    bool
    is_hetalt() const
    {
        unsigned print_gt(smod.max_gt);
        const uint8_t a0(DIGT::get_allele(print_gt,0));
        const uint8_t a1(DIGT::get_allele(print_gt,1));
        return ((a0!=a1) && (dgt.ref_gt != a0) && (dgt.ref_gt != a1));
    }

    bool
    is_nonref() const override
    {
        return (smod.max_gt != dgt.ref_gt);
    }

    bool
    is_print_unknowngt() const
    {
        return (smod.is_unknown || (!smod.is_used_covered));
    }

    bool
    is_deletion() const
    {
        return ((!smod.is_unknown) && smod.is_used_covered && (!smod.is_zero_ploidy) && (is_nonref()));
    }

    bool
    is_qual() const
    {
        return ((!smod.is_unknown) && smod.is_used_covered && (!smod.is_zero_ploidy) && (is_nonref()));
    }

    std::string phased_ref, phased_alt, phased_AD;
    diploid_genotype dgt;
    double hapscore = 0;
    double MQ = 0;				 // RMS of mapping qualities

    //only meaningful for het calls
    double ReadPosRankSum = 0;  // Uses Mann-Whitney Rank Sum Test for the distance from the end of the read containing an alternate allele.
    double BaseQRankSum = 0;    // Uses Mann-Whitney Rank Sum Test for BQs (ref bases vs alternate alleles)
    double MQRankSum = 0;       // Uses Mann-Whitney Rank Sum Test for MQs (ref bases vs alternate alleles)
    double avgBaseQ = 0;
    double rawPos = 0;
    unsigned mapq_zero = 0;     // The number of spanning reads that do not pass the command-line mapq test

    digt_call_info smod;
};

std::ostream& operator<<(std::ostream& os,const digt_site_info& si);

struct continuous_site_info : public site_info
{
    continuous_site_info(const pos_t init_pos,
                         const char init_ref,
                         const snp_pos_info& good_pi,
                         const int used_allele_count_min_qscore,
                         const double min_het_vf,
                         const bool is_forced_output = false) : site_info(init_pos,
                                                                              init_ref,
                                                                              good_pi,
                                                                              used_allele_count_min_qscore,
                                                                              is_forced_output)
        , _min_het_vf(min_het_vf)
    {
    }

    bool is_snp() const override
    {
        return _is_snp;
    }
    void set_filter (VCF_FILTERS::index_t filter) override
    {
        for (auto& call : calls) call.set_filter(filter);
    }
    bool is_nonref() const override
    {
        auto ref_id = base_to_id(ref);
        return calls.end() !=
               std::find_if(calls.begin(), calls.end(),
                            [&](const continuous_site_call& call)
        {
            return call._base != ref_id;
        });
    }

    bool _is_snp = false;


    const char* get_gt(const continuous_site_call& call) const
    {
        if (call._base == base_to_id(ref))
            return "0/0";
        else if (call.variant_frequency() >= (1 -_min_het_vf))
            return "1/1";
        else if (call.variant_frequency() < _min_het_vf)
            return "0/0"; // STAR-66 - desired behavior
        else
            return "0/1";
    }

    std::vector<continuous_site_call> calls;
private:
    double _min_het_vf;
};

struct continuous_indel_info : public indel_info
{
    explicit continuous_indel_info(const pos_t init_pos)
        : indel_info(init_pos)
    {
    }

    void set_filter (VCF_FILTERS::index_t filter) override
    {
        for (auto& call : calls) call.set_filter(filter);
    }


    bool is_indel() const override
    {
        for (auto& call : calls)
        {
            if (call._alleleDepth > 0)
                return true;
        }
        return false;
    }

    bool is_forced_output() const override
    {
        for (auto& call : calls)
            if (call._id.is_forced_output)
                return true;
        return false;
    }

    pos_t end() const override
    {
        pos_t result = 0;
        for (auto& x : calls)
            result = std::max(result, x._ik.right_pos());
        return result;
    }

    const char* get_gt() const
    {
        if (is_het)
            return "0/1";
        else
            return "1/1";
    }

    std::vector<continuous_indel_call> calls;
    bool is_het=false;


};



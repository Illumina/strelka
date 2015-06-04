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



struct shared_modifiers
{
    shared_modifiers()
    {
        clear();
    }

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
    }

    int gqx;
    int gq;
    unsigned max_gt;

    std::bitset<VCF_FILTERS::SIZE> filters;
};


std::ostream& operator<<(std::ostream& os,const shared_modifiers& shmod);


struct indel_modifiers : public shared_modifiers
{
    indel_modifiers()
    {
        clear();
    }

    void
    clear()
    {
        shared_modifiers::clear();
        ploidy.clear();
        cigar.clear();
        Qscore = -1;
    }

    ALIGNPATH::path_t cigar;

    /// represent site ploidy over the reference span of the overlapping indel set in the event of overlap:
    std::vector<unsigned> ploidy;

    // The empirically calibrated quality-score of an indel, if -1 no q-score has been reported
    int Qscore = -1;
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

struct site_modifiers : public shared_modifiers
{
    site_modifiers()
    {
        clear();
    }

    void
    clear()
    {
        shared_modifiers::clear();
        is_unknown=true;
        is_covered=false;
        is_used_covered=false;
        is_zero_ploidy=false;
        is_block=false;
        is_phased_region=false;
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
    bool is_block;
    bool is_phased_region;

    MODIFIED_SITE_GT::index_t modified_gt;

    // The empirically calibrated quality-score of the site, if -1 not q-score has been reported
    int Qscore = -1;
};


std::ostream& operator<<(std::ostream& os,const site_modifiers& smod);



struct indel_info
{
    indel_info(const pos_t init_pos,
         const indel_key& init_ik,
         const starling_diploid_indel_core& init_dindel,
         const starling_indel_report_info& init_iri,
         const starling_indel_sample_report_info& init_isri)
    {
        pos=(init_pos);
        ik=(init_ik);
        dindel=(init_dindel);
        _iri.push_back(init_iri);
        _isri.push_back(init_isri);
        _imod.emplace_back();
    }

    void add_overlap(const reference_contig_segment& ref, indel_info& overlap);

    const char*
    get_gt() const
    {
        if (this->is_hetalt())
        {
            return "1/2";
        }
        else if (dindel.is_haploid())
        {
            using namespace STAR_DIINDEL;

            switch (imod().max_gt)
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
        return STAR_DIINDEL::get_gt_label(imod().max_gt);
    }

    bool
    is_hetalt() const
    {
        return (_imod.size() > 1);
    }

    bool
    is_het() const
    {
        return (static_cast<int>(imod().max_gt)>1);
    }

    // the site ploidy within the indel at offset x
    unsigned
    get_ploidy(const unsigned offset) const
    {
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
            assert(offset<imod().ploidy.size());
            return imod().ploidy[offset];
        }
        return 2;
    }
    void set_hap_cigar(
        const unsigned lead=1,
        const unsigned trail=0);

    std::map<std::string, double>
    get_indel_qscore_features(const double chrom_depth) const;

    pos_t pos;
    indel_key ik;
    starling_diploid_indel_core dindel;
    std::vector<starling_indel_report_info> _iri;
    std::vector<starling_indel_sample_report_info> _isri;
    std::vector<indel_modifiers> _imod;

    indel_modifiers& imod() { return _imod.front(); }
    const indel_modifiers& imod(unsigned idx=0) const { return _imod[idx]; }

    //starling_indel_report_info& iri() { return _iri.front(); }
    const starling_indel_report_info& iri(unsigned idx = 0) const { return _iri[idx]; }

    //starling_indel_sample_report_info& isri() { return _isri.front(); }
    const starling_indel_sample_report_info& isri(unsigned idx=0) const { return _isri[idx]; }



};


//Data structure defining parameters for a single site to be used for writing in gvcf_aggregator
struct site_info
{
    site_info()
    {
        std::fill(known_counts.begin(),known_counts.end(),0);
    }

    void
    init(const pos_t init_pos,
         const char init_ref,
         const snp_pos_info& good_pi,
         const int used_allele_count_min_qscore)
    {
        pos=(init_pos);
        ref=(init_ref);
        good_pi.get_known_counts(known_counts,used_allele_count_min_qscore);

        phased_ref.clear();
        phased_alt.clear();
        phased_AD.clear();

        Unphasable = false;

        smod.clear();
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
            unsigned print_gt(smod.max_gt);
            if (smod.is_block)
            {
                print_gt = dgt.ref_gt;
            }
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
        return DIGT::is_het(print_gt) && ref != DIGT::label(print_gt)[0] && ref != DIGT::label(print_gt)[1];
    }

    bool
    is_nonref() const
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
        return ((!smod.is_block) && (!smod.is_unknown) && smod.is_used_covered && (!smod.is_zero_ploidy) && (is_nonref()));
    }

    bool
    is_qual() const
    {
        return ((!smod.is_block) && (!smod.is_unknown) && smod.is_used_covered && (!smod.is_zero_ploidy) && (is_nonref()));
    }

    pos_t pos = 0;
    char ref = 'N';
    std::string phased_ref, phased_alt, phased_AD;
    unsigned n_used_calls = 0;
    unsigned n_unused_calls = 0;
    std::array<unsigned,N_BASE> known_counts;
    diploid_genotype dgt;
    unsigned hpol = 0;
    double hapscore = 0;
    double MQ = 0;				 // RMS of mapping qualities

    //only meaningful for het calls
    double ReadPosRankSum = 0;  // Uses Mann-Whitney Rank Sum Test for the distance from the end of the read containing an alternate allele.
    double BaseQRankSum = 0;    // Uses Mann-Whitney Rank Sum Test for BQs (ref bases vs alternate alleles)
    double MQRankSum = 0;       // Uses Mann-Whitney Rank Sum Test for MQs (ref bases vs alternate alleles)
    double avgBaseQ = 0;
    double rawPos = 0;
    unsigned mapq_zero = 0;     // The number of spanning reads that do not pass the command-line mapq test
    bool Unphasable = false;        // Set to true if the site should never be included in a phasing block

    site_modifiers smod;
};

std::ostream& operator<<(std::ostream& os,const site_info& si);

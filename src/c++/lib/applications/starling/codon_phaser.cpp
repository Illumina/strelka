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
/*
 * Codon_phaser.cpp
 *
 *  Created on: Sep 10, 2013
 *  Author: Morten Kallberg
 */

#include "codon_phaser.hh"

#include <array>
#include <sstream>
#include <vector>

//#define DEBUG_CODON


#ifdef DEBUG_CODON
#include "blt_util/log.hh"
#endif


// Add a SNP site to the phasing buffer
bool
Codon_phaser::
add_site(const site_info& si)
{
    _buffer.push_back(si);
#ifdef DEBUG_CODON
    log_os << __FUNCTION__ << ": input si " << si << "\n";
#endif

    // case: extending block with het call, update block_end position
    if (is_phasable_site(si))
    {
        if (! is_in_block())
        {
            block_start = si.pos;
#ifdef DEBUG_CODON
            log_os << __FUNCTION__ << ": phasable & starting block @ " << block_start << "\n";
#endif
        }
#ifdef DEBUG_CODON
        else
        {
            log_os << __FUNCTION__ << ": phasable & continuing block to " << si.pos << "\n";
        }
#endif
        block_end = si.pos;
        het_count++;
        return false;
    }

    //case: we get a record that is explicitly set a unphaseable
    if (si.Unphasable)
    {
#ifdef DEBUG_CODON
        log_os << __FUNCTION__ << ": I shouldn't phase this record " << si << "\n";
#endif
        return true;
    }

    // case: extending block with none-het call based on the phasing range
    if (is_in_block() && (si.pos-block_end+1)<this->opt.phasing_window)
    {
#ifdef DEBUG_CODON
        log_os << __FUNCTION__ << ": notphasable & continuing block to " << si.pos << "\n";
#endif
        return false;
    }

    // case: setup the phased record and write out
    if (het_count>1)
    {
//        log_os << "!!!het count " << het_count << "\n";
#ifdef DEBUG_CODON
//        this->write_out_buffer();
#endif
        make_record();
    }
    return true;
}

void
Codon_phaser::construct_reference()
{
    this->reference = "";
    for (unsigned i=0; i<this->_buffer.size()-(this->opt.phasing_window-1); i++)
        this->reference += _buffer.at(i).ref;
}

void
Codon_phaser::
create_phased_record()
{
    // sanity check do we have all record we expect or are there un-accounted no-calls
    if (this->get_block_length()>this->_buffer.size())
        return;

    if (total_reads<10)
    {
        // some initial minimum conditions, look for at least 10 spanning reads support
        // set flag on records saying too little evidence to phase
        //        for (int i=0;i<this->buffer.size();i++)
        //            if(this->buffer.at(i).is_het())
        //                this->buffer.at(i).smod.;
        return;
    }

    assert(_buffer.size()>0);

    //Decide if we accept the novel alleles, very hacky criteria for now
    //at least 80% of the reads have to support a diploid model
    //TODO unphased genotype corresponds to as many phased genotypes as there are permutations of
    //the observed alleles. Thus, for a given unphased genotyping G1, . . . ,Gn,
    //we need to to calculate probability of sampling a particular unphased genotype combinations
    //given a set of allele frequencies...
    typedef std::pair<std::string, int> allele_count_t;
    //static constexpr unsigned alleleCount(2);
    std::array<allele_count_t, 2> max_alleles = {allele_count_t("N",0),allele_count_t("N",0)};
    for (auto& obs : observations)
    {
        if (obs.second>max_alleles[0].second)
        {
            max_alleles[1]  = max_alleles[0];
            max_alleles[0]  = obs;
        }
        if (obs.second>max_alleles[1].second && max_alleles[0].first!=obs.first)
        {
            max_alleles[1] = obs;
        }
    }

#ifdef DEBUG_CODON
    log_os << "max_1 " << max_alleles[0].first << "=" << max_alleles[0].second << "\n";
    log_os << "max_2 " << max_alleles[1].first << "=" << max_alleles[1].second << "\n";
#endif

    // some ad hoc metrics to measure consistency with diploid model
    const int allele_sum = max_alleles[0].second + max_alleles[1].second;
    const float max_allele_frac = (1.0*allele_sum)/this->total_reads;
    const float relative_allele_frac  = 1.0*max_alleles[1].second/max_alleles[0].second;

#ifdef DEBUG_CODON
    log_os << "max alleles sum " << allele_sum << "\n";
    log_os << "max alleles frac " << max_allele_frac << "\n";
    log_os << "relative_allele_frac " << relative_allele_frac << "\n";
#endif

    bool phasing_inconsistent(false);
    if (max_allele_frac<0.8)
    {
        // non-diploid?
        phasing_inconsistent = true;
#ifdef DEBUG_CODON
        log_os << __FUNCTION__ << "; non-diploid\n";
#endif
    }

    if (relative_allele_frac<0.5)
    {
        // allele imbalance?
        phasing_inconsistent = true;
#ifdef DEBUG_CODON
        log_os << __FUNCTION__ << "; allele imbalance\n";
#endif
    }

    // sanity check that we have one het on each end of the block:
    {
        assert(! max_alleles[0].first.empty());
        assert(max_alleles[0].first.size() == max_alleles[1].first.size());
        if (max_alleles[0].first.front() == max_alleles[1].first.front()) phasing_inconsistent=true;
        if (max_alleles[0].first.back() == max_alleles[1].first.back()) phasing_inconsistent=true;
#ifdef DEBUG_CODON
        log_os << __FUNCTION__ << "; non-sane\n";
#endif
    }

    if (phasing_inconsistent)
    {
        for (auto& val : _buffer)
        {
            if (! is_phasable_site(val)) continue;
            val.smod.set_filter(VCF_FILTERS::PhasingConflict);
        }
        return;
    }

    // we have a phased record, modify site buffer to reflect the changes
    site_info& base = (this->_buffer.at(0));

    base.phased_ref = this->reference;
    const bool is_ref(max_alleles[0].first==this->reference || max_alleles[1].first==this->reference);

    base.smod.is_block = false;
    base.smod.is_unknown = false;
    base.smod.max_gt = 4;
    base.dgt.ref_gt  = 0; // hacking  the gt method to 0/1
    if (!is_ref) base.dgt.ref_gt = 2; // hacking the gt method to 1/2

    std::stringstream AD,alt;
    AD << this->observations[this->reference];
    for (const auto& val : max_alleles)
    {
        if (val.first==reference) continue;
        if (! alt.str().empty()) alt << ',';
        alt << val.first;
        AD << ',' << val.second;
    }

#ifdef DEBUG_CODON
    log_os << "buffer size " << _buffer.size() << "\n";
    log_os << "block length " << get_block_length() << "\n";
#endif

    // set GQ and GQX
    static const int maxInt(std::numeric_limits<int>::max());
    int min_gq(maxInt), min_qual(maxInt), min_qscore(maxInt);
    for (unsigned i(0); i<this->get_block_length(); i++)
    {
        const site_info& si(_buffer.at(i));
        if (! is_phasable_site(si)) continue;
        min_gq = std::min(si.smod.gq,min_gq);
        min_qual = std::min(si.dgt.genome.snp_qphred,min_qual);
        min_qscore = std::min(si.smod.Qscore,min_qscore);
    }

    // set various quality fields conservatively
    base.smod.gq                = min_gq;
    base.dgt.genome.snp_qphred  = min_qual;
    base.smod.gqx               = std::min(min_gq,min_qual);
    base.smod.Qscore            = min_qscore;

    base.phased_alt = alt.str();
    base.phased_AD  = AD.str();

    base.smod.is_phased_region = true;
    const int reads_ignored = (this->total_reads-allele_sum);
    base.n_used_calls = this->total_reads - reads_ignored;
    base.n_unused_calls = this->total_reads_unused + reads_ignored; // second term mark all alleles that we didnt use as unused reads

    // case we want to report the phased record clean out buffer and push on the phased record only
#ifdef DEBUG_CODON
    log_os << "buffer size " << _buffer.size() << "\n";
    this->write_out_buffer();
#endif

    // erase buffer sites which are now combined in the buffer[0] record
    _buffer.erase(_buffer.begin()+1,_buffer.begin()+this->get_block_length());
}

// makes the phased VCF record from the buffered sites list
void
Codon_phaser::
make_record()
{
    this->construct_reference();
    this->collect_pileup_evidence();
    this->create_phased_record();
}



void
Codon_phaser::
collect_pileup_evidence()
{
    // build quick pileup index over phase range:
    std::vector<const snp_pos_info*> spi;
    for (int blockPos(block_start); blockPos<=block_end; ++blockPos)
    {
        const snp_pos_info& pi(bc_buff.get_pos(blockPos));
        spi.push_back(&pi);
    }

    const unsigned blockWidth(spi.size());
    std::vector<int> callOffset(blockWidth,0);

    // analyze as virtual reads -- to do so treat the first pileup column as a privilaged refernece point:
    const snp_pos_info& beginPi(*spi[0]);
    const unsigned n_calls(beginPi.calls.size());
    for (unsigned callIndex(0); callIndex<n_calls; ++callIndex)
    {
        std::string sub_str;
        bool isPass(true);
        for (unsigned blockIndex(0);blockIndex<blockWidth;blockIndex++)
        {
            const snp_pos_info& pi(*spi[blockIndex]);
            while (true)
            {
                const base_call& bc(pi.calls[callIndex+callOffset[blockIndex]]);
                if (bc.isLastBaseCallFromMatchSeg)
                {
                    for (unsigned blockIndex2(blockIndex+1);blockIndex2<blockWidth;++blockIndex2)
                    {
                        callOffset[blockIndex2]--;
                    }
                }
                if (bc.isFirstBaseCallFromMatchSeg && (blockIndex != 0))
                {
                    for (unsigned blockIndex2(blockIndex);blockIndex2<blockWidth;++blockIndex2)
                    {
                        callOffset[blockIndex2]++;
                    }
                }
                else
                {
                    break;
                }
            }
            const base_call& bc(pi.calls[callIndex+callOffset[blockIndex]]);
            if (bc.is_call_filter)
            {
                isPass=false;
            }
            if (bc.isLastBaseCallFromMatchSeg && ((blockIndex+1) < blockWidth))
            {
                isPass=false;
                break;
            }
            sub_str += id_to_base(bc.base_id);
        }
#ifdef DEBUG_CODON
        log_os << __FUNCTION__ << ": callOffset:";
        for (const auto off : callOffset)
        {
            log_os << "\t" << off;
        }
        log_os << "\n";
        log_os << __FUNCTION__ << ": callIndex, isPass, substr: " << callIndex << " " << isPass << " " << sub_str << "\n";
#endif

        if (! isPass) continue;

        observations[sub_str]++;
        total_reads++;
    }

    total_reads_unused=n_calls-total_reads;
}



void
Codon_phaser::clear()
{
    _buffer.clear();
    observations.clear();
    block_start = -1;
    block_end   = -1;
    het_count                   = 0;
    total_reads                 = 0;
    total_reads_unused          = 0;
    reference                   = "";
}



void
Codon_phaser::write_out_buffer() const
{
    for (const auto& val : _buffer)
    {
        log_os << val << " ref " << val.ref << "\n";
    }
}

void
Codon_phaser::write_out_alleles() const
{
    for (const auto& val : observations)
    {
        if (val.first == this->reference) continue;
        log_os << val.first << "=" << val.second << " ";
    }
    log_os << "\n";
}

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
/*
 *
 *  Created on: Sep 10, 2013
 *  Author: Morten Kallberg
 */

#include "codon_phaser.hh"

#include <array>
#include <functional>
#include <vector>

//#define DEBUG_CODON


#ifdef DEBUG_CODON
#include "blt_util/log.hh"
#endif


void Codon_phaser::process(std::unique_ptr<GermlineSiteLocusInfo> site)
{
    std::unique_ptr<GermlineDiploidSiteLocusInfo> si(downcast<GermlineDiploidSiteLocusInfo>(std::move(site)));
    if (opt.do_codon_phasing && (is_phasable_site(si) || is_in_block()))
    {
        auto emptyBuffer = add_site(std::move(si));
        if ((!is_in_block()) || emptyBuffer)
        {
            output_buffer();
        }
        return;
    }
    _sink->process(std::move(si));
}

void Codon_phaser::flush_impl()
{
    if (opt.do_codon_phasing)
    {
        if (het_count>1)
        {
            make_record();
        }
        output_buffer();
    }
}

// the Codon phaser can't work across indels, so flush any in-progress phasing
void Codon_phaser::process(std::unique_ptr<GermlineIndelLocusInfo> ii)
{
    if (opt.do_codon_phasing && is_in_block())
    {
        if (het_count>1)
        {
            make_record();
        }
        output_buffer();
    }
    _sink->process(std::move(ii));
}

void Codon_phaser::output_buffer()
{
    for (auto& si : _buffer)
    {
        _sink->process(std::move(si));
    }
    clear();
}


// Add a SNP site to the phasing buffer
bool
Codon_phaser::add_site(std::unique_ptr<GermlineDiploidSiteLocusInfo> si)
{
#ifdef DEBUG_CODON
    log_os << __FUNCTION__ << ": input si " << *si << "\n";
#endif

    // case: extending block with het call, update block_end position
    if (is_phasable_site(si))
    {
        if (! is_in_block())
        {
            block_start = si->pos;
#ifdef DEBUG_CODON
            log_os << __FUNCTION__ << ": phasable & starting block @ " << block_start << "\n";
#endif
        }
#ifdef DEBUG_CODON
        else
        {
            log_os << __FUNCTION__ << ": phasable & continuing block to " << si->pos << "\n";
        }
#endif
        block_end = si->pos;
        het_count++;
        _buffer.push_back(std::move(si));
        return false;
    }

    //case: we get a record that is explicitly set a unphaseable
    if (si->Unphasable)
    {
#ifdef DEBUG_CODON
        log_os << __FUNCTION__ << ": I shouldn't phase this record " << *si << "\n";
#endif
        _buffer.push_back(std::move(si));

        return true;
    }

    // case: extending block with none-het call based on the phasing range
    if (is_in_block() && (si->pos-block_end+1)<this->opt.phasing_window)
    {
#ifdef DEBUG_CODON
        log_os << __FUNCTION__ << ": notphasable & continuing block to " << si->pos << "\n";
#endif
        _buffer.push_back(std::move(si));

        return false;
    }
    _buffer.push_back(std::move(si));


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
    reference = "";
    auto start = _buffer.front()->pos;
    auto end = _buffer.back()->pos;
    auto len = std::min(get_block_length(), (unsigned)(end-start+1));

    ref.get_substring(start, len, reference);

#ifdef DEBUG_CODON
    log_os << __FUNCTION__ << ": reference " << reference << "\n";
#endif
}

static bool is_genotype_represented(int genotype, unsigned offset, const std::string& allele1, const std::string& allele2)
{
    const char* gt = DIGT::label(genotype);
    for (int i = 0; i < 2; i++)
    {
        char base = gt[i];
        if ( (allele1.length() <= offset || allele1[offset] != base) &&
             (allele2.length() <= offset || allele2[offset] != base))
        {
            return false;
        }
    }
    return true;
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
        for (auto& si : _buffer)
        {
            if (si->is_het())
                si->smod.is_phasing_insufficient_depth = true;
        }
        return;
    }

    assert(_buffer.size()>0);



    //Decide if we accept the novel alleles, very hacky criteria for now
    //at least 80% of the reads have to support a diploid model
    //TODO unphased genotype corresponds to as many phased genotypes as there are permutations of
    //the observed alleles. Thus, for a given unphased genotyping G1, . . . ,Gn,
    //we need to to calculate probability of sampling a particular unphased genotype combinations
    //given a set of allele frequencies...
    // TODO: PASSing SNPs are getting filtered during phasing for allele balance issues. Looks like it would benefit from a
    // model-based approach.
    // < chr8  70364286    .   G   A   32  PASS    SNVSB=-2.5;SNVHPOL=4    GT:GQ:GQX:DP:DPF:AD 0/1:65:21:15:5:12,3
    // < chr8  70364287    .   A   G   11  PASS    SNVSB=-2.5;SNVHPOL=2    GT:GQ:GQX:DP:DPF:AD 0/1:44:18:16:3:14,2
    //---
    //> chr8  70364286    .   G   A   32  PhasingConflict SNVSB=-2.5;SNVHPOL=4    GT:GQ:GQX:DP:DPF:AD 0/1:65:21:15:5:12,3
    //> chr8  70364287    .   A   G   11  PhasingConflict SNVSB=-2.5;SNVHPOL=2    GT:GQ:GQX:DP:DPF:AD 0/1:44:18:16:3:14,2
    //static constexpr unsigned alleleCount(2);
    std::array<std::pair<std::string,allele_observations>, 2> max_alleles;
    for (const auto& obs : observations)
    {
#ifdef DEBUG_CODON
        log_os << "obs:" << obs.first << "(" << obs.second.count() << ")\n";
#endif

        if (obs.second.count()>max_alleles[0].second.count())
        {
            max_alleles[1] = max_alleles[0];
            max_alleles[0] = obs;
        }
        if (obs.second.count()>max_alleles[1].second.count() && max_alleles[0].first!=obs.first)
        {
            max_alleles[1] = obs;
        }
    }

#ifdef DEBUG_CODON
    log_os << "max_1 " << max_alleles[0].first << "=" << max_alleles[0].second.count() << "\n";
    log_os << "max_2 " << max_alleles[1].first << "=" << max_alleles[1].second.count() << "\n";
#endif

    // some ad hoc metrics to measure consistency with diploid model
    const int allele_sum = max_alleles[0].second.count() + max_alleles[1].second.count();
    const float max_allele_frac = (static_cast<float>(allele_sum))/this->total_reads;
    const float relative_allele_frac = static_cast<float>(max_alleles[1].second.count())/max_alleles[0].second.count();

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
    if (! phasing_inconsistent)
    {
        if (max_alleles[0].second.count() == 0) phasing_inconsistent=true;
        if (max_alleles[1].second.count() == 0) phasing_inconsistent=true;

        if (! phasing_inconsistent)
        {
            assert(! max_alleles[0].first.empty());
            assert(max_alleles[0].first.size() == max_alleles[1].first.size());
            if (max_alleles[0].first.front() == max_alleles[1].first.front()) phasing_inconsistent=true;
            if (max_alleles[0].first.back() == max_alleles[1].first.back()) phasing_inconsistent=true;
        }

#ifdef DEBUG_CODON
        if (phasing_inconsistent) log_os << __FUNCTION__ << "; non-sane\n";
#endif
    }

    std::stringstream AD, ADF, ADR, alt;
    if (!phasing_inconsistent)
    {
        AD << observations[reference].count();
        ADF << observations[reference].fwd;
        ADR << observations[reference].rev;
        for (const auto& val : max_alleles)
        {
            if (val.first==reference) continue;
            if (! alt.str().empty()) alt << ',';
            alt << val.first;
            AD << ',' << val.second.count();
            ADF << ',' << val.second.fwd;
            ADR << ',' << val.second.rev;
        }
    }

#ifdef DEBUG_CODON
    log_os << "buffer size " << _buffer.size() << "\n";
    log_os << "block length " << get_block_length() << "\n";
#endif

    // set GQ and GQX
    static const int maxInt(std::numeric_limits<int>::max());
    int min_qual(maxInt), min_EVS(maxInt);
    std::vector<unsigned> pls;
    unsigned ref_gt(0);
    unsigned max_gt(0);
    bool is_min_gq_idx0(false);
    bool is_min_gq_idx1(false);
    unsigned min_gq_idx0(0);
    unsigned min_gq_idx1(0);
    if (!phasing_inconsistent)
    {
        for (unsigned i(0); i<this->get_block_length(); i++)
        {
            const auto& si(_buffer.at(i));
            if (! is_phasable_site(si)) continue;
            // It is possible for a weak hetalt call and a strong het call to be phased into a het
            // call, particularly at a site that is triallelic with one ref allele. In that case, it is
            // important to skip those calls' contribution to the resultant statistics. This logic causes
            // us to skip over these variants.
            if (!is_genotype_represented(si->smod.max_gt, unsigned(si->pos - _buffer[0]->pos),
                                         max_alleles[0].first, max_alleles[1].first))
            {
#ifdef DEBUG_CODON
                log_os << "Variant at offset i lacks support in phased alleles - failing phasing\n"
#endif
                       // TODO: This is a blunt instrument. Refine this solution to select the minimum set of
                       // ALTs that support the buffered variants
                       phasing_inconsistent = true;
                continue;
            }

            if ((! is_min_gq_idx0) || (si->smod.gq < _buffer.at(min_gq_idx0)->smod.gq))
            {
                min_gq_idx1 = min_gq_idx0;
                if (is_min_gq_idx0) is_min_gq_idx1 = true;
                min_gq_idx0 = i;
                is_min_gq_idx0 = true;
            }
            if ((i != min_gq_idx0) && ( (! is_min_gq_idx1) || (si->smod.gq < _buffer.at(min_gq_idx1)->smod.gq)))
            {
                min_gq_idx1 = i;
                is_min_gq_idx1 = true;
            }
            min_qual = std::min(si->dgt.genome.snp_qphred,min_qual);
            min_EVS = std::min(si->empiricalVariantScore,min_EVS);
        }
    }
    if (!is_min_gq_idx0)
    {
#ifdef DEBUG_CODON
        log_os << __FUNCTION__ << "; no variants left\n";
#endif
        phasing_inconsistent = true;
    }


    if (phasing_inconsistent)
    {
        for (auto& val : _buffer)
        {
            if (! is_phasable_site(val)) continue;
            val->filters.set(GERMLINE_VARIANT_VCF_FILTERS::PhasingConflict);
        }
        return;
    }

    const bool is_ref(max_alleles[0].first==this->reference || max_alleles[1].first==this->reference);

    int min_gq(maxInt);
    {
        const auto& minsi0(*(_buffer.at(min_gq_idx0)));
        min_gq = minsi0.smod.gq;
        max_gt = minsi0.smod.max_gt;
        pls = minsi0.dgt.phredLoghood;
        ref_gt = minsi0.dgt.ref_gt;

        if (! is_ref && is_min_gq_idx1)
        {
            //
            // hetalt case
            //
#ifdef DEBUG_CODON
            log_os << "min0 " << min_gq_idx0 << "\n";
            log_os << "min1 " << min_gq_idx1 << "\n";
#endif

            // create fake ref_gt value for hetalt case:
            const uint8_t ax(DIGT::get_allele(max_gt,0));
            const uint8_t ay(DIGT::get_allele(max_gt,1));

            // phase a0/a1 to match max_allele order;
            const bool is_swap(max_alleles[1].first[min_gq_idx0] == id_to_base(ax));
            const uint8_t a0(is_swap ? ay : ax);
            const uint8_t a1(is_swap ? ax : ay);

            ref_gt = 0;
            for (; true; ref_gt++)
            {
                assert(ref_gt < N_BASE);
                if ((ref_gt != a0) && (ref_gt != a1)) break;
            }

            // hetalt site:
            const auto& minsi1(*(_buffer.at(min_gq_idx1)));
            const uint8_t bx(DIGT::get_allele(minsi1.smod.max_gt,0));
            const uint8_t by(DIGT::get_allele(minsi1.smod.max_gt,1));

            // phase b0/b1 to match max_allele order;
            const bool is_swap2(max_alleles[1].first[min_gq_idx1] == id_to_base(bx));
            const uint8_t b0(is_swap2 ? by : bx);
            const uint8_t b1(is_swap2 ? bx : by);

#ifdef DEBUG_CODON
            log_os << "ref/a0/a1/b0/b1 " << ref_gt << " " << (int)a0 << " " << (int)a1 << " " << (int)b0 << " " << (int)b1 << "\n";
            log_os << "ref0/ref1 " << minsi0.dgt.ref_gt << " " << minsi1.dgt.ref_gt << "\n";
#endif

            // construct new fake approximated PL distribution:
            const auto& pls0(minsi0.dgt.phredLoghood);
            const auto& pls1(minsi1.dgt.phredLoghood);

            pls[ref_gt] = pls0[minsi0.dgt.ref_gt] + pls1[minsi1.dgt.ref_gt];  // 0/0
            pls[DIGT::get_gt_with_alleles(ref_gt,a0)] = pls0[DIGT::get_gt_with_alleles(minsi0.dgt.ref_gt,a0)] + pls1[DIGT::get_gt_with_alleles(minsi1.dgt.ref_gt,b0)];  // 0/1
            pls[DIGT::get_gt_with_alleles(a0,a0)] = pls0[DIGT::get_gt_with_alleles(a0,a0)] + pls1[DIGT::get_gt_with_alleles(b0,b0)];  // 1/1
            pls[DIGT::get_gt_with_alleles(ref_gt,a1)] = pls0[DIGT::get_gt_with_alleles(minsi0.dgt.ref_gt,a1)] + pls1[DIGT::get_gt_with_alleles(minsi1.dgt.ref_gt,b1)];  // 0/2
            pls[max_gt] = 0;  // 1/2
            pls[DIGT::get_gt_with_alleles(a1,a1)] = pls0[DIGT::get_gt_with_alleles(a1,a1)] + pls1[DIGT::get_gt_with_alleles(b1,b1)];  // 2/2
        }
    }

    // we have a phased record, modify site buffer to reflect the changes
    auto& base = this->_buffer.at(0);

    base->phased_ref = this->reference;
    base->smod.is_unknown = false;
    base->smod.max_gt = max_gt;
    base->dgt.ref_gt = ref_gt;

    // set various quality fields conservatively
    base->smod.gq                = min_gq;
    base->dgt.genome.snp_qphred  = min_qual;
    base->dgt.phredLoghood       = pls;
    base->smod.gqx               = std::min(min_gq,min_qual);
    base->empiricalVariantScore  = min_EVS;

    base->phased_alt = alt.str();
    base->phased_AD  = AD.str();
    base->phased_ADF  = ADF.str();
    base->phased_ADR  = ADR.str();

    base->smod.is_phased_region = true;
    const int reads_ignored = (this->total_reads-allele_sum);
    base->n_used_calls = this->total_reads - reads_ignored;
    base->n_unused_calls = this->total_reads_unused + reads_ignored; // second term mark all alleles that we didnt use as unused reads

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

    /// traces individual read fragments from the pileup structure within the phasing block range
    ///
    /// Low detail summary:
    /// If translating aligned reads into pileup columns is thought of as a sort of matrix
    /// transpose, this function is trying to invert the transposition back to a (partial)
    /// read. A naive pileup structure would not support this, but starling pileup information
    /// has been supplemented with a few extra bits that allows this reconstruction.
    ///
    /// Given a start offset within the phasing block (startBlockIndex), return a tuple composed of:
    ///     1. A bool indicating whether the read is good (ie. complete to the end of the
    ///        phasing block and all basecalls pass filter) This return val is really only
    ///        useful for startBlockIndex=0
    ///     2. The reconstructed read fragment
    ///
    /// External calls should always provide startBlockIndex=0, other values of this argument
    /// are used by internal recursive calls to the function. If the return bool indicates a
    /// good read, this function effectively returns a single allele count for the phaser
    ///
    /// Important implementation detail: this function mutates the external vector 'callOffset'
    /// over successive calls to progressively jump to the offset of the next read
    ///
    /// isFirstBaseCallFromMatchSeg is set for any basecall starting a continuous matching sequence in one read
    /// isLastBaseCallFromMatchSeg is set for any basecall ending a continuous matching sequence in one read
    // can't use auto for the return value here b/c of recursion:
    std::function<std::pair<bool,std::string>(const unsigned)> tracePartialRead =
        [&](const unsigned startBlockIndex)
    {
        bool isPass(true);
        std::string readFrag;
        for (unsigned blockIndex(startBlockIndex); blockIndex<blockWidth; ++blockIndex)
        {
            const snp_pos_info& pi(*spi[blockIndex]);
            while (true)
            {
                const base_call& bc(pi.calls.at(callOffset[blockIndex]));
                if ((! bc.isFirstBaseCallFromMatchSeg) || (blockIndex == startBlockIndex)) break;
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wuninitialized"
                tracePartialRead(blockIndex);
#pragma clang diagnostic pop
            }
            const base_call& bc(pi.calls.at(callOffset[blockIndex]));
            // this represents the 'ordinary' advance of callOffset for a non-interrupted read segment:
            callOffset[blockIndex]++;
            if (bc.isLastBaseCallFromMatchSeg && ((blockIndex+1) < blockWidth))
            {
                isPass=false;
                break;
            }
            if (bc.is_call_filter) isPass=false;
            readFrag += id_to_base(bc.base_id);
        }
        return std::make_pair(isPass,readFrag);
    };

    // analyze as virtual reads -- to do so treat the first pileup column as a privileged reference point:
    const snp_pos_info& beginPi(*spi[0]);
    const unsigned n_calls(beginPi.calls.size());
    for (unsigned callIndex(0); callIndex<n_calls; ++callIndex)
    {
        const bool is_fwd_strand(beginPi.calls[callIndex].is_fwd_strand);
        std::pair<bool,std::string> result = tracePartialRead(0);
#ifdef DEBUG_CODON
        log_os << __FUNCTION__ << ": callOffset:";
        for (const auto off : callOffset)
        {
            log_os << "\t" << off;
        }
        log_os << "\n";
        log_os << __FUNCTION__ << ": callIndex, isPass, substr: " << callIndex << " " << result.first << " " << result.second << "\n";
#endif

        if (! result.first) continue;

        if (is_fwd_strand)
        {
            observations[result.second].fwd++;
        }
        else
        {
            observations[result.second].rev++;
        }
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
Codon_phaser::
write_out_buffer(std::ostream& os) const
{
    for (const auto& val : _buffer)
    {
        os << *val << " ref " << val->ref << "\n";
    }
}

void
Codon_phaser::
write_out_alleles(std::ostream& os) const
{
    for (const auto& val : observations)
    {
        if (val.first == this->reference) continue;
        os << val.first << "=" << val.second.count() << " ";
    }
    os << "\n";
}

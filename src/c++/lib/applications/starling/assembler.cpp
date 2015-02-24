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

#include "assembler.hh"

#include <array>
#include <sstream>
#include <vector>

//#define DEBUG_ASSEMBLE


#ifdef DEBUG_ASSEMBLE
#include "blt_util/log.hh"
#endif


// Add a SNP site to the phasing buffer
bool
assembler::
add_site(const site_info& si)
{
    _buffer.push_back(si);

    // case: extending block with variant call
    if (si.is_nonref())
    {
        if (! is_in_block())
            block_start = si.pos;
        block_end = si.pos;
        het_count ++;
#ifdef DEBUG_ASSEMBLE
        log_os << "starting block @ " << (this->block_start+1) << " with " << si << "\n";
#endif
        return false;
    }

    //case: we get a record that is explicitly set to not be assembled
    if (si.Unphasable)
    {
#ifdef DEBUG_ASSEMBLE
        log_os << "I shouldn't assemble this record " << si << "\n";
#endif
        return true;
    }

    // case: extending block with none-het call based on the phasing range
    if (is_in_block() && (si.pos-block_end+1)<this->opt.phasing_window)
    {
#ifdef DEBUG_ASSEMBLE
        log_os << "Extending block with @ " << (this->block_start+1) << " with " << si << "\n";
#endif
        return false;
    }

    // case: setup the assembled records
    if (het_count>1)
    {
//        log_os << "!!!record count " << het_count << "\n";
#ifdef DEBUG_ASSEMBLE
//        this->write_out_buffer();
#endif
        make_record();
    }
    return true;
}

void
assembler::construct_reference()
{
    this->reference = "";
    for (unsigned i=0; i<this->_buffer.size()-(this->opt.phasing_window-1); i++)
        this->reference += _buffer.at(i).ref;
}

void
assembler::create_contig_record()
{
    // sanity check do we have all record we expect or are there un-accounted no-calls
    if (this->get_block_length()>this->_buffer.size())
        return;

    if (this->total_reads<10)
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

#ifdef DEBUG_ASSEMBLE
    log_os << "max_1 " << max_alleles[0].first << "=" << max_alleles[0].second << "\n";
    log_os << "max_2 " << max_alleles[1].first << "=" << max_alleles[1].second << "\n";
#endif

    // some ad hoc metrics to measure consistency with diploid model
    const int allele_sum = max_alleles[0].second + max_alleles[1].second;
    const float max_allele_frac = (1.0*allele_sum)/this->total_reads;
    const float relative_allele_frac  = 1.0*max_alleles[1].second/max_alleles[0].second;

#ifdef DEBUG_ASSEMBLE
    log_os << "max alleles sum " << allele_sum << "\n";
    log_os << "max alleles frac " << max_allele_frac << "\n";
    log_os << "relative_allele_frac " << relative_allele_frac << "\n";
#endif

    bool phasing_inconsistent(false);
    if (max_allele_frac<0.8)
    {
        // non-diploid?
        phasing_inconsistent = true;
    }

    if (relative_allele_frac<0.5)
    {
        // allele imbalance?
        phasing_inconsistent = true;
    }

    if (phasing_inconsistent)
    {
        for (auto& val : _buffer)
        {
            if (! val.is_het()) continue;
            val.smod.set_filter(VCF_FILTERS::PhasingConflict); // switch to custom filter for assembly conlfict.
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

#ifdef DEBUG_ASSEMBLE
    log_os << "buffer size " << _buffer.size() << "\n";
    log_os << "block length " << get_block_length() << "\n";
#endif

    // set GQ and GQX
    int min_gq(INT_MAX), min_qual(INT_MAX), min_qscore(INT_MAX);
    for (unsigned i(0); i<this->get_block_length(); i++)
    {
        const site_info& si(_buffer.at(i));
        if (! si.is_het()) continue;
        min_gq = std::min(si.smod.gq,min_gq);
        min_qual = std::min(si.dgt.genome.snp_qphred,min_qual);
        min_qscore = std::min(si.Qscore,min_qscore);
    }

    // set various quality fields conservatively
    base.smod.gq                = min_gq;
    base.dgt.genome.snp_qphred  = min_qual;
    base.smod.gqx               = std::min(min_gq,min_qual);
    base.Qscore                 = min_qscore;

    base.phased_alt = alt.str();
    base.phased_AD  = AD.str();

    base.smod.is_phased_region = true;
    int reads_ignored = (this->total_reads-allele_sum);
    base.n_used_calls = this->total_reads - reads_ignored;
    base.n_unused_calls = this->total_reads_unused + reads_ignored; // second term mark all alleles that we didnt use as unused reads

    // case we want to report the phased record clean out buffer and push on the phased record only
#ifdef DEBUG_ASSEMBLE
    log_os << "buffer size " << _buffer.size() << "\n";
    this->write_out_buffer();
#endif

    // erase buffer sites which are now combined in the buffer[0] record
    _buffer.erase(_buffer.begin()+1,_buffer.begin()+this->get_block_length());
}

// makes the phased VCF record from the buffered sites list
void
assembler::make_record()
{
    this->construct_reference();
    this->collect_read_evidence();
    this->create_contig_record();
}

void
assembler::
collect_read_segment_evidence(
    const read_segment& rseg)
{
    //check if we are covering the block range
    if (static_cast<unsigned>(abs(rseg.buffer_pos-this->block_end))>rseg.read_size())
        return;

    const bam_seq bseq(rseg.get_bam_read());

    // read quality checks
    if (static_cast<int>(rseg.map_qual())<this->opt.min_single_align_score || rseg.is_invalid_realignment || !rseg.is_valid())
    {
//                this->total_reads_unused++; // do not count filtered reads in DPF
        return;
    }

    const int sub_start((this->block_start-rseg.buffer_pos));
    const int sub_end((this->block_end-rseg.buffer_pos));
#ifdef DEBUG_ASSEMBLE
    int pad(0); // add context to extracted alleles for debugging
    sub_start -= pad;
    sub_end += pad;
#endif
    using namespace BAM_BASE;

    if (sub_start>=0 && sub_end<max_read_len)
    {
        bool do_include(true);
        //instead make bit array counting the number of times het pos co-occur
        std::string sub_str("");
        for (int t=sub_start; t<(sub_end+1); t++) //pull out substring of read
        {
            if (bseq.get_char(t)=='N'|| static_cast<int>(rseg.qual()[t]<this->opt.min_qscore))
            {
                do_include = false; // do qual check of individual bases here, kick out entire read if we dont meet cut-off
                break;
            }
            sub_str+= bseq.get_char(t); //TODO use more efficient data structure here
        }
#ifdef DEBUG_ASSEMBLE
//                    log_os << "substart " << sub_start << "\n";
//                    log_os << "subend " << sub_end << "\n";
//                    log_os << "substr " << sub_str << "\n";
//                    log_os << "read key " << rseg.key() << "\n";
//                    log_os << "read pos " << rseg.buffer_pos << "\n";
//                    log_os << "do_include " << do_include << "\n";
//                    log_os << "read seq " << bseq << "\n\n";
#endif
        if (do_include)
        {
            this->observations[sub_str]++;
            total_reads++;
        }
        else
            this->total_reads_unused++;
    }
}



void
assembler::
collect_read_evidence()
{
    int buffer_pos = (block_start-this->max_read_len);
    const int buffer_end = (block_start);
    // extract evidence for all reads that span the entire phasing range
    for (; buffer_pos<buffer_end; buffer_pos++)
    {
        read_segment_iter ri(read_buffer.get_pos_read_segment_iter(buffer_pos));
        read_segment_iter::ret_val r;
        while (true)
        {
            r=ri.get_ptr();
            if (nullptr==r.first) break;
            collect_read_segment_evidence(r.first->get_segment(r.second));
            ri.next();
        }
    }

#ifdef DEBUG_ASSEMBLE
    log_os << "max read len " << max_read_len << "\n";
#endif

#ifdef DEBUG_ASSEMBLE
    log_os << "total reads " << total_reads << "\n";
    log_os << "total reads unused " << total_reads_unused << "\n";
#endif
}

void
assembler::clear()
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
assembler::write_out_buffer() const
{
    for (const auto& val : _buffer)
    {
        log_os << val << " ref " << val.ref << "\n";
    }
}

void
assembler::write_out_alleles() const
{
    for (const auto& val : observations)
    {
        if (val.first == this->reference) continue;
        log_os << val.first << "=" << val.second << " ";
    }
    log_os << "\n";
}

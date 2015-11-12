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

#include "alignment_util.hh"
#include "starling_read_align_score_indels.hh"

#include "blt_util/log.hh"
#include "starling_common/indel_util.hh"

#include <algorithm>
#include <iostream>


//#define DEBUG_ALIGN


typedef std::pair<indel_key,std::pair<bool,indel_key> > indel_status_t;
typedef std::pair<double,const candidate_alignment*> align_info_t;
typedef std::map<indel_status_t,align_info_t> iks_map_t;



// check whether the current path has a better score than that stored
// and if so add entry to iks map
//
static
void
check_and_update_iks(iks_map_t& iks_map,
                     const indel_key& ik_call,
                     const bool is_indel_present,
                     const indel_key& ik_present,
                     const double path_lnp,
                     const candidate_alignment* cal_ptr)
{
    const indel_status_t mkey(std::make_pair(ik_call,std::make_pair(is_indel_present,ik_present)));

    // check to see if the path score is better than what we already have:
    const iks_map_t::const_iterator j(iks_map.find(mkey));
    if ((j!=iks_map.end()) && ((j->second.first) >= path_lnp)) return;
    iks_map[mkey]=std::make_pair(path_lnp,cal_ptr);
}



typedef std::map<indel_key,indel_set_t> overlap_map_t;



static
void
overlap_map_tick(overlap_map_t& omap,
                 const indel_key& ik1,
                 const indel_key& ik2)
{
    indel_set_t& os(omap[ik1]);
    if (os.find(ik2) == os.end()) os.insert(ik2);
}



static
bool
is_interfering_indel(const indel_set_t& current_indels,
                     const indel_key& new_indel)
{
    if (current_indels.count(new_indel) != 0) return false;

    for (const indel_key& ik : current_indels)
    {
        if (is_indel_conflict(ik,new_indel)) return true;
    }
    return false;
}



// by how many bases does the alignment cross each indel breakpoint?
//
// Note that this function only works correctly for indels that are
// present in the alignment -- otherwise there are multiple possible
// realignments of the read
//
// Allow for the possibility of an upstream oligo to be implicitely
// present for the read, acting to 'anchor' the beginning of the read
// with the equivalent length match state even though this does not
// appear in the read's CIGAR string.
//
static
std::pair<int,int>
get_alignment_indel_bp_overlap(const unsigned upstream_oligo_size,
                               const alignment& al,
                               const indel_key& ik)
{
    using namespace ALIGNPATH;

    // all read positions are expressed in unclipped read coordinates!
    //
    pos_t read_head_pos(0);
    pos_t ref_head_pos(al.pos);

    // mark the left and right indel breakpoints in read position
    // coordinates:
    bool is_left_read_pos(false);
    pos_t left_read_pos(0);
    bool is_right_read_pos(false);
    pos_t right_read_pos(0);

    for (const path_segment& ps : al.path)
    {
        pos_t next_read_head_pos(read_head_pos);
        pos_t next_ref_head_pos(ref_head_pos);

        if       (is_segment_align_match(ps.type))
        {
            next_read_head_pos += ps.length;
            next_ref_head_pos += ps.length;
        }
        else if (ps.type == INSERT)
        {
            next_read_head_pos += ps.length;
        }
        else if (ps.type == DELETE)
        {
            next_ref_head_pos += ps.length;
        }
        else if ((ps.type == SOFT_CLIP) || (ps.type == HARD_CLIP))
        {
            // do nothing... this function operates on unclipped read
            // coordinates instead of true read coordinates
            //
        }
        else
        {
            // TODO deal with other CIGAR types....
            assert(0);
        }

        if ((! is_left_read_pos) && (ik.pos<=(next_ref_head_pos)))
        {
            left_read_pos=read_head_pos+(ik.pos-ref_head_pos);
            is_left_read_pos=true;
        }
        if ((! is_right_read_pos) && (ik.right_pos()<(next_ref_head_pos)))
        {
            right_read_pos=read_head_pos+(ik.right_pos()-ref_head_pos);
            is_right_read_pos=true;
        }

        read_head_pos=next_read_head_pos;
        ref_head_pos=next_ref_head_pos;
    }

    // extension values account for oligo anchors at the end of the
    // read not represented in the CIGAR string:
    int left_extension(0);
    int right_extension(0);
    if (al.is_fwd_strand)
    {
        if ( left_read_pos > 0 ) left_extension = upstream_oligo_size;
    }
    else
    {
        if ( (read_head_pos-right_read_pos)  > 0 )right_extension = upstream_oligo_size;
    }

    int left_overlap(0);
    if (is_left_read_pos)
    {
        left_overlap=std::min(left_read_pos+left_extension,(read_head_pos-left_read_pos));
        left_overlap=std::max(0,left_overlap);
    }

    int right_overlap(0);
    if (is_right_read_pos)
    {
        right_overlap=std::min(right_read_pos,(read_head_pos-right_read_pos)+right_extension);
        right_overlap=std::max(0,right_overlap);
    }

    return std::make_pair(left_overlap,right_overlap);
}



typedef std::set<std::pair<indel_key,indel_key> > indel_pair_set;



// given two alignemnts that have (nearly) the same alignment score,
// quickly approximate whether they could be expressing the same set
// of indels:
//
static
bool
is_equiv_candidate(const candidate_alignment& cal1,
                   const candidate_alignment& cal2,
                   const unsigned max_indel_size,
                   indel_pair_set& equiv_keys)
{
    equiv_keys.clear();

    indel_set_t is1,is2;
    get_alignment_indels(cal1,max_indel_size,is1);
    get_alignment_indels(cal2,max_indel_size,is2);

    const unsigned s1(is1.size());
    const unsigned s2(is2.size());

    if (s1 != s2) return false;

    indel_set_t::const_iterator it1(is1.begin()), it1_end(is1.end());
    indel_set_t::const_iterator it2(is2.begin()); //, it2_end(is2.end());

    for (; it1!=it1_end; ++it1,++it2)
    {
        if (*it1==*it2) continue;
        if (it1->type != it2->type) return false;
        if (it1->length != it2->length) return false;
        equiv_keys.insert(std::make_pair(*it1,*it2));
    }
    return true;
}



// written for equiv indel resolution (i.e. assuming ik1 and ik2 are
// of the same type and length)
//
static
bool
is_first_indel_dominant(
    const indel_synchronizer& isync,
    const indel_key& ik1,
    const indel_key& ik2)
{
    const bool ic1(isync.is_candidate_indel(ik1));
    const bool ic2(isync.is_candidate_indel(ik2));

    if (ic2 && (! ic1)) return false;
    if (ic2==ic1)
    {
        return (ik1.pos<=ik2.pos);
    }
    return true;
}



// use the most likely alignment for each indel state for every indel
// in indel_status_map to generate data needed in indel calling:
//
void
score_indels(const starling_base_options& opt,
             const starling_base_deriv_options& /*dopt*/,
             const starling_sample_options& sample_opt,
             const read_segment& rseg,
             indel_synchronizer& isync,
             const std::set<candidate_alignment>& cal_set,
             const bool is_incomplete_search,
             const std::vector<double>& cal_set_path_lnp,
             double max_path_lnp,
             const candidate_alignment* max_cal_ptr)
{
    static const bool is_safe_mode(true);

    // (1) score candidate alignments -- already done before calling this function


    // (2) store the highest scoring alignment with each indel present
    // and absent.
    //
    // Note there are up to two 'indel absent' states. The first is
    // when an interfering indel is not present, and the second is
    // when an interfering indel is present. The second case will be
    // less common.
    //
    // We must control for the case where the indel is absent but
    // replaced by an equivilent indel. In this case we must not only
    // determine that an equivilent exists, but establish which case
    // is dominant so that the indel is only reported once. Dominance
    // rules should follow the same pattern established for
    // realignment above.
    //
    // The indel normalization here is intended to be a final fallback
    // corrective measure to prevent indels from disappearing at the
    // call stage. It is expected (and preferred) that as much
    // upstream indel normalization will take place as is practical.
    //

    // (2alpha)
    //
    // go through all scored alignments and search for cases with an
    // equal score (or within the smooth threshold if this is turned
    // on). Determine if smooth score matches are actually equivilent
    // -- in which case remove the non-dominant alignment and mark the
    // equivilent indel as not participating in this read's indel call
    // search:
    //
    const unsigned cal_set_size(cal_set.size());

    // Nonnorm indels contains any indels found to be equivilent to
    // another indel (ie. not normalized). These indels are not used
    // in subsequent calculation so as not to penalize the normalized
    // indel form (otherwise all reads which support the indel would
    // appear to ambiguously support two different indels).
    //
    indel_set_t nonnorm_indels;

    std::vector<bool> cal_set_exclude(cal_set_size,false);

    // This procedure implements a heuristic late-stage indel
    // normalization, set is_slip_norm=false to disable:
    //
    static const bool is_slip_norm(true);
    if (is_slip_norm)
    {

        const double equiv_lnp_range( opt.is_smoothed_alignments ?
                                      opt.smoothed_lnp_range : 0. );

        // go through alignment x alignments in best->worst score
        // order -- to do so start out with the sort order:
        std::vector<std::pair<double,unsigned> > sorted_path_lnp;
        std::vector<const candidate_alignment*> cal_ptr_vec;
        {
            std::set<candidate_alignment>::const_iterator si(cal_set.begin());
            for (unsigned i(0); i<cal_set_size; ++i,++si)
            {
                sorted_path_lnp.push_back(std::make_pair(cal_set_path_lnp[i],i));
                cal_ptr_vec.push_back(&(*si));
            }
            std::sort(sorted_path_lnp.rbegin(),sorted_path_lnp.rend());
        }

        // this starts out as a copy of cal_set_path_lnp, but path_lnp
        // gets reset to the highest lnp an alignment has eliminated
        // to prevent transitive shenanigans:
        std::vector<double> smooth_path_lnp(cal_set_path_lnp);

        bool is_any_excluded(false);
        indel_pair_set ips;

        for (unsigned i(0); i<cal_set_size; ++i)
        {
            const unsigned sorti(sorted_path_lnp[i].second);
            if (cal_set_exclude[sorti]) continue;

            for (unsigned j(i+1); j<cal_set_size; ++j)
            {
                const unsigned sortj(sorted_path_lnp[j].second);
                if (cal_set_exclude[sortj]) continue;

                // determine if this is an equiv_indel:
                if (smooth_path_lnp[sortj]+equiv_lnp_range < smooth_path_lnp[sorti]) break;
                const bool is_equiv(is_equiv_candidate(*(cal_ptr_vec[sorti]),
                                                       *(cal_ptr_vec[sortj]),
                                                       opt.max_indel_size,
                                                       ips));

                if (! is_equiv) continue;

#ifdef DEBUG_ALIGN
                log_os << "COWSLIP: sorti,lnpi,sortj,lnpj: "
                       << sorti << " " << smooth_path_lnp[sorti] << " "
                       << sortj << " " << smooth_path_lnp[sortj] << "\n";
                log_os << "COWSLIP cali: " << *(cal_ptr_vec[sorti]) << "\n";
                log_os << "COWSLIP calj: " << *(cal_ptr_vec[sortj]) << "\n";
#endif

                // In theory, this should be an assertion of
                // !ips.empty(), but it seems that two alignments with
                // the same score and the same set of indels, but
                // (i.e. different start positions) could somehow
                // occur.
                if (ips.empty()) continue;

                indel_pair_set::const_iterator ip(ips.begin()),ip_end(ips.end());
                // for each pair, determine which indel is dominant,
                // if there are multiple equiv indels, we remove the
                // alignment containing the first non-dominant equiv
                // indel only
                bool is_sorti_removed(false);
                bool is_removed(false);
                for (; ip!=ip_end; ++ip)
                {
                    const bool is1(is_first_indel_dominant(isync,ip->first,ip->second));

#ifdef DEBUG_ALIGN
                    log_os << "COWSLIP: indel1: " << ip->first << "\n";
                    log_os << "COWSLIP: indel2: " << ip->second << "\n";
                    log_os << "COWSLIP: is_indel1_dominant?: " << is1 << "\n";
#endif
                    if (is1)
                    {
                        nonnorm_indels.insert(ip->second);
#ifdef DEBUG_ALIGN
                        log_os << "COWSLIP: marking 2 nonnorm: " << ip->second << "\n";
#endif
                        if (! is_removed)
                        {
                            cal_set_exclude[sortj] = true;
                            is_any_excluded = true;
                            smooth_path_lnp[sorti] = std::max(smooth_path_lnp[sorti],smooth_path_lnp[sortj]);
#ifdef DEBUG_ALIGN
                            log_os << "COWSLIP: excluding sortj: " << sortj << "\n";
#endif
                        }
                    }
                    else
                    {
                        nonnorm_indels.insert(ip->first);
#ifdef DEBUG_ALIGN
                        log_os << "COWSLIP: marking 1 nonnorm: " << ip->first << "\n";
#endif
                        if (! is_removed)
                        {
                            cal_set_exclude[sorti] = true;
                            is_any_excluded = true;
                            smooth_path_lnp[sortj] = std::max(smooth_path_lnp[sorti],smooth_path_lnp[sortj]);
#ifdef DEBUG_ALIGN
                            log_os << "COWSLIP: excluding sorti: " << sorti << "\n";
#endif
                            is_sorti_removed=true;
                        }
                    }
                    is_removed=true;
                }
                if (is_sorti_removed) break;
            }
        }

        // recalc new max path in case it was excluded:
        //
        // note this loop is designed to take advantage of the high->low path_lnp sort
        //
        if (is_any_excluded)
        {
            for (unsigned i(0); i<cal_set_size; ++i)
            {
                const unsigned sorti(sorted_path_lnp[i].second);
                if (cal_set_exclude[sorti]) continue;
#ifdef DEBUG_ALIGN
                log_os << "COWSLIP: reseting max_path_lnp/max_cal from: " << max_path_lnp << " " << *(max_cal_ptr) << "\n";
#endif
                max_path_lnp=cal_set_path_lnp[sorti];
                max_cal_ptr=cal_ptr_vec[sorti];
#ifdef DEBUG_ALIGN
                log_os << "COWSLIP: reseting max_path_lnp/max_cal to: " << max_path_lnp << " " << *(max_cal_ptr) << "\n";
#endif
                break;
            }
        }
    }


    // (2a) get an initial set of candidate indels which can be scored
    //      from max path. note that we can't do full
    //      breakpoint-testing yet, so a few more indels may be
    //      eliminated from this set:
    //
    // The criteria for indel calling are that:
    // 1) the indel is not private to the read
    // 2) the read overlaps at least one indel breakpoint by at least option:min_read_bp_flank bases
    // 3) the indel is not in the nonnorm set
    // 4) the read meets whatever mapping quality criteria have been setup for snp-calling
    //
    // graveyard:
    // //?? 2) the most likely alignment path intersects the indel    //
    // //?? 4) the most likely read alignment does not contain an interfering indel
    //
    // TODO deal with mapping issues somehow??
    //
    const bool is_tier1_read(rseg.is_treated_as_tier1_mapping());

    const candidate_alignment& max_cal(*max_cal_ptr);

    indel_set_t max_cal_eval_indels;
    indel_set_t max_cal_indels;
    indel_status_map_t is_max_cal_eval_indels_interfere;
    {
        const known_pos_range max_pr(get_soft_clip_alignment_range(max_cal.al));
        const known_pos_range strict_max_pr(get_strict_alignment_range(max_cal.al));

#ifdef DEBUG_ALIGN
        log_os << "VARMIT: starting max_cal_eval_indel search\n"
               << "VARMIT: max_cal.al: " << max_cal.al << "\n"
               << "VARMIT: max_pr: " << max_pr << "\n"
               << "VARMIT: strict_max_pr: " << strict_max_pr << "\n";
#endif
        get_alignment_indels(max_cal,opt.max_indel_size,max_cal_indels);
#ifdef DEBUG_ALIGN
        log_os << "VARMIT max_path extracted indels:\n";
        dump_indel_set(max_cal_indels,log_os);
#endif
        indel_buffer& ibuff(isync.ibuff());
        const std::pair<iiter,iiter> ipair(ibuff.pos_range_iter(max_pr.begin_pos,max_pr.end_pos));
        for (iiter i(ipair.first); i!=ipair.second; ++i)
        {
            const indel_key& ik(i->first);
            indel_data& id(get_indel_data(i));

#ifdef DEBUG_ALIGN
            log_os << "VARMIT: max path eval indel candidate: " << ik;
#endif

            if (! isync.is_candidate_indel(ik,id)) continue;

#ifdef DEBUG_ALIGN
            log_os << "VARMIT: max path indel is candidate\n";
#endif

            const bool is_indel_interfering(is_interfering_indel(max_cal_indels,ik));

#ifdef DEBUG_ALIGN
            log_os << "VARMIT: max path indel is non-interfering\n";
#endif

            const bool is_indel_present(max_cal_indels.count(ik)!=0);
#ifdef DEBUG_ALIGN
            log_os << "VARMIT: indel present? " << is_indel_present << "\n";
#endif

            if (is_indel_present)
            {
                const std::pair<int,int> both_bpo(get_alignment_indel_bp_overlap(opt.upstream_oligo_size,max_cal.al,ik));
                const int bpo(std::max(both_bpo.first,both_bpo.second));
#ifdef DEBUG_ALIGN
                log_os << "VARMIT: indel bp_overlap " << bpo << "\n";
#endif
                if (bpo < sample_opt.min_read_bp_flank)
                {
                    if (bpo>0)
                    {
                        if (is_tier1_read) id.suboverlap_tier1_read_ids.insert(rseg.id());
                        else              id.suboverlap_tier2_read_ids.insert(rseg.id());
                    }
                    continue;
                }
            }
            else     // check that the most likely alignment intersects this indel:
            {
#ifdef DEBUG_ALIGN
                log_os << "VARMIT: indel intersects max_path? "
                       <<  is_range_intersect_indel_breakpoints(strict_max_pr,ik) << "\n";
#endif
                if (! is_range_intersect_indel_breakpoints(strict_max_pr,ik)) continue;
            }

            // all checks passed... indel will be evaluated for indel calling:
            //
            max_cal_eval_indels.insert(ik);

            is_max_cal_eval_indels_interfere[ik] = is_indel_interfering;
        }
    }

    // go through eval indel set and map which indels conflict with
    // each other
    //
    // for now we calc and store this info in a comically inefficient
    // manner
    //
    overlap_map_t indel_overlap_map;
    {
        indel_set_t::const_iterator i(max_cal_eval_indels.begin());
        const indel_set_t::const_iterator i_end(max_cal_eval_indels.end());
        for (; i!=i_end; ++i)
        {
            indel_set_t::const_iterator j(i);
            ++j;
            for (; j!=i_end; ++j)
            {
                if (is_indel_conflict(*i,*j))
                {
                    overlap_map_tick(indel_overlap_map,*i,*j);
                    overlap_map_tick(indel_overlap_map,*j,*i);
                }
            }
        }
    }



    // (2b) Create index of highest scoring alignment set for each
    //      evaluation indel (found in step 2a). The set of alignments
    //      include all mutually consistent sets of indels which
    //      interfere with the evalutation indel (including the eval
    //      indel itself) 99% of the time this should just be the
    //      indel and ref paths.  99% of the rest of the time this
    //      will just be the paths of the indel, ref, and one
    //      alternate indel allele.
    //

    //
    // First part of iks_map_t key is the called indel -- second part
    // of the key is the indel which is toggled on which interferes
    // with the called indel. When called indel is toggled on the
    // first and second keys are the same. Ptr is NULL for no-indel
    // (reference) allele case.
    //
    // Effect of this change is now we store maximum scoring alignment
    // for each interfering indel allele (not just the called indel
    // and reference). I don't realistically envision that more than
    // one indel will be required very often -- that's the case where we
    // would want a full-blown candidate haplotype model anyway.
    //
    //
    iks_map_t iks_max_path_lnp;
    {
        {
            // as a simple acceleration to the full max scoring
            // alignment path search below -- we go through all indel
            // states present in the global maximum scoring alignment,
            // taking advantage of our knowledge that this will be the
            // highest scoring path already:
            //
            const align_info_t max_info(std::make_pair(max_path_lnp,&max_cal));
            for (const indel_key& eval_ik : max_cal_eval_indels)
            {
                const bool is_indel_present(max_cal_indels.count(eval_ik)!=0);

                if (is_indel_present)
                {
                    iks_max_path_lnp[std::make_pair(eval_ik,std::make_pair(is_indel_present,eval_ik))] = max_info;

                    // mark this as an alternate indel score for interfering indels:
                    for (const indel_key& overlap_ik : indel_overlap_map[eval_ik])
                    {
                        iks_max_path_lnp[std::make_pair(overlap_ik,std::make_pair(is_indel_present,eval_ik))] = max_info;
                    }
                }
                else
                {
                    // check that this indel does not interfere with the max-set:
                    if (! is_max_cal_eval_indels_interfere[eval_ik])
                    {
                        iks_max_path_lnp[std::make_pair(eval_ik,std::make_pair(is_indel_present,eval_ik))] = max_info;
                    }
                }
            }
        }

        std::set<candidate_alignment>::const_iterator si(cal_set.begin()),si_end(cal_set.end());
        for (unsigned c(0); si!=si_end; ++si,++c)
        {
            const candidate_alignment& ical(*si);
            const bool is_max_cal(&ical == &max_cal);

            if (cal_set_exclude[c])
            {
                assert(! is_max_cal);
                continue;
            }
            if (is_max_cal) continue;

            const double path_lnp(cal_set_path_lnp[c]);

            indel_set_t ical_indels;
            get_alignment_indels(ical,opt.max_indel_size,ical_indels);

            for (const indel_key& eval_ik : max_cal_eval_indels)
            {
                const bool is_indel_present(ical_indels.count(eval_ik)!=0);
                if (is_indel_present)
                {
                    //const indel_status_t mkey(std::make_pair(eval_ik,std::make_pair(is_indel_present,eval_ik)));
                    check_and_update_iks(iks_max_path_lnp,eval_ik,is_indel_present,eval_ik,path_lnp,&ical);

                    // mark this as an alternate indel score for interfering indels:
                    for (const indel_key& overlap_ik : indel_overlap_map[eval_ik])
                    {
                        check_and_update_iks(iks_max_path_lnp,overlap_ik,is_indel_present,eval_ik,path_lnp,&ical);
                    }
                }
                else
                {
                    // if indel is not present, we must determine
                    // whether an interfering indel is present in the
                    // alignment to score this alignment in the right
                    // category:
                    //
                    bool is_interference(false);
                    for (const indel_key& overlap_ik : indel_overlap_map[eval_ik])
                    {
                        if (ical_indels.count(overlap_ik)!=0)
                        {
                            is_interference=true;
                            break;
                        }
                    }

                    if (! is_interference)
                    {
                        check_and_update_iks(iks_max_path_lnp,eval_ik,is_indel_present,eval_ik,path_lnp,&ical);
                    }
                }
            }
        }
    }


    // (3) assemble indel calling scores for each indel into the indel buffer:
    //
    //
    const unsigned read_length(rseg.read_size());

    // for indel caller calculation we always need the "rest-of-genome"
    // alignment score, which is calculated from the number of non-ambiguous sites:
    uint16_t nsite(0);
    {
        const bam_seq bseq(rseg.get_bam_read());
        for (unsigned i(0); i<read_length; ++i)
        {
            if (bseq.get_code(i) == BAM_BASE::ANY) continue;
            nsite++;
        }
    }


    // Note that for sanity/consistency, the breakpoint overlap
    // filtration criteria is applied to the called indel
    // only. Alignments of the same read to alternate indel alleles
    // are not tested for overlap size
    //
    {
        indel_buffer& ibuff(isync.ibuff());
        for (const indel_key& eval_ik : max_cal_eval_indels)
        {
            indel_data* id_ptr(ibuff.get_indel_data_ptr(eval_ik));
            assert(NULL != id_ptr);

            // we test for presence of the indel on the highest
            // scoring alignment because breakpoint overlap has
            // already been tested there, allowing us to skip this
            // step:
            const bool is_indel_present_on_max_path(max_cal_indels.count(eval_ik)!=0);

#ifdef DEBUG_ALIGN
            log_os << "VARMIT: final indel scan: " << eval_ik;
            log_os << "VARMIT: is_present_on_max_path?: " << is_indel_present_on_max_path << "\n";
#endif

            // sanity check that these exist in iks (unless search is incomplete):
            // 1) indel present
            // 2) indel absent
            // 3) conflicting indels present
            //

            // 1) indel present:
            double indel_path_lnp(max_path_lnp);
            if (! is_indel_present_on_max_path)
            {
                const iks_map_t::iterator j(iks_max_path_lnp.find(std::make_pair(eval_ik,std::make_pair(true,eval_ik))));
                const bool is_found(j!=iks_max_path_lnp.end());
                if (! is_found)
                {
                    if (is_incomplete_search) continue;
                    // TODO -- get more precise information on exactly when we expect an indel which overlaps a max_cal indel to not be found, for
                    // now we have to give a pass on all cases:
                    if (is_max_cal_eval_indels_interfere[eval_ik]) continue;

                    // for the nonnorm cases, we've already eliminated at least some of the alignments which contain them
                    // and in certain circumstances there won't be alternates available
                    if (nonnorm_indels.count(eval_ik)!=0) continue;

                    if (is_safe_mode)
                    {
                        log_os << "WARNING: ";
                    }
                    else
                    {
                        log_os << "ERROR: ";
                    }

                    log_os << "failed to find expected alignment for indel: " << eval_ik
                           << "\twhile evaluating read_segment:\n" << rseg << "\n";

                    if (is_safe_mode)
                    {
                        continue;
                    }
                    else
                    {
                        exit(EXIT_FAILURE);
                    }
                }

                const candidate_alignment& alt_cal(*(j->second.second));
                const std::pair<int,int> both_bpo(get_alignment_indel_bp_overlap(opt.upstream_oligo_size,alt_cal.al,eval_ik));
                const int bpo(std::max(both_bpo.first,both_bpo.second));

#ifdef DEBUG_ALIGN
                log_os << "VARMIT: called_indel_bp_overlap " << bpo << "\n";
#endif
                if (bpo < sample_opt.min_read_bp_flank)
                {
                    if (bpo>0)
                    {
                        if (is_tier1_read) id_ptr->suboverlap_tier1_read_ids.insert(rseg.id());
                        else              id_ptr->suboverlap_tier2_read_ids.insert(rseg.id());
                    }
                    continue;
                }

                indel_path_lnp=j->second.first;
            }

#ifdef DEBUG_ALIGN
            log_os << "VARMIT: called_indel_path_lnp " << indel_path_lnp << "\n";
#endif

            // 2) indel absent w/o interference:
            double ref_path_lnp(0);
            {
                const iks_map_t::iterator j(iks_max_path_lnp.find(std::make_pair(eval_ik,std::make_pair(false,eval_ik))));
                const bool is_found(j!=iks_max_path_lnp.end());
                if (is_incomplete_search && (! is_found)) continue;
                if (! is_found)
                {
                    if (is_safe_mode)
                    {
                        log_os << "WARNING: ";
                    }
                    else
                    {
                        log_os << "ERROR: ";
                    }

                    log_os << "failed to find reference alignment while evaluating read_segment:\n" << rseg << "\n";

                    if (is_safe_mode)
                    {
                        continue;
                    }
                    else
                    {
                        exit(EXIT_FAILURE);
                    }
                }

                ref_path_lnp=j->second.first;

#ifdef DEBUG_ALIGN
                log_os << "VARMIT: ref_path_lnp " << ref_path_lnp << "\n";
#endif
            }

            // assemble basic score data:
            //
            int readpos = max_cal.al.is_fwd_strand ? ((signed) (eval_ik.pos - max_cal.al.pos)) :
                          ((signed) (read_length - eval_ik.pos + max_cal.al.pos - 1));

#ifdef DEBUG_ALIGN
            //            correct for strandedness
            log_os << "indelpos: " << eval_ik.pos << " readpos: + " << rseg.genome_align().pos << " rp: " << readpos << "\n";
#endif
            read_path_scores rps(ref_path_lnp,indel_path_lnp,nsite,read_length,is_tier1_read,max_cal.al.is_fwd_strand,
                                 (int16_t) readpos);

            // start adding alternate indel alleles, if present:

            for (const indel_key& overlap_ik : indel_overlap_map[eval_ik])
            {
                const iks_map_t::iterator j(iks_max_path_lnp.find(std::make_pair(eval_ik,std::make_pair(true,overlap_ik))));

#ifdef DEBUG_ALIGN
                log_os << "VARMIT: alternate_indel " << overlap_ik;
#endif

                // TODO consider a way that the basic stores
                // could still be used even if we can't get all
                // alternate paths?
                //
                const bool is_found(j!=iks_max_path_lnp.end());

                //if(is_incomplete_search and (not is_found)) continue;
                //assert(is_found);

                // for now, this case has to be allowed, because it is possible
                // for no alignment not including the called indel to be in range
                // of the alternate indel
                //
                // for now we continue when the alt case is not found.
                //
                // TODO tighten the re-alignment procedure so that we know when the
                // alt case should be missing -- possibly use ref_path instead of continuing
                // in this case
                //
                if (is_found)
                {
                    rps.insert_alt(overlap_ik,j->second.first);
                    //                        rps.alt_indel[*k] = j->second.first;
                }
                else
                {
                    // TODO look at these cases in more detail before allowing them!
                    continue;
                    //rp.alt_indel[*k] = ref_path_lnp;
                }

#ifdef DEBUG_ALIGN
                //log_os << "VARMIT: alternate_indel_path_found: " << is_found << "\n";
                log_os << "VARMIT: alternate_indel_path_lnp " << j->second.first << "\n";
#endif
            }

            id_ptr->read_path_lnp[rseg.id()] = rps;
            id_ptr->n_mapq++;
            id_ptr->cumm_mapq += rseg.map_qual();
            if(rseg.map_qual() == 0)
            {
                id_ptr->n_mapq0++;
            }

#ifdef DEBUG_ALIGN
            log_os << "VARMIT: modified indel data: " << *(id_ptr);
#endif
        }
    }
}


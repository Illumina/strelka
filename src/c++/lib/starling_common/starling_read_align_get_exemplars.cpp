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

/// \file
///
/// \author Chris Saunders
///

#include "blt_util/log.hh"
#include "starling_common/alignment_util.hh"
#include "starling_common/starling_read_align_get_exemplars.hh"

#include <cassert>

#include <iostream>


//#define DEBUG_ALIGN



// return true if alignments are compatible -- the only allowed
// difference is edge soft clipping, if a match exists in the new
// alignment where the exemplar is clipped, the match will be
// transfered
//
// For now we assume both alignments have been "matchified" -- will
// need better soft-clip policy later.
//
// Without more information we have to assume edge inserts are private
// to the discovery alignment. This means that two otherwise
// compatible discovery alignments separated by an edge insertions are
// incompatible.  A next iteration of these routines could come up
// with a better private indel handling policy. Ideally we would want
// some step to predetermine what was private in the discovery
// alignment beforehand? But even non-private indels need to be
// preserved in the D.A. so that we no which indel alleles to evaluate
// for indel scoring -- given that the indel eval algorithm remains
// constrained to two alleles only.
//
static
bool
check_and_adjust_exemplar(const alignment& al1,
                          alignment& al2)
{
    return (al1==al2);
}






// test new alignment against the exemplar set to determine if it
// is compatible with one that already exists:
//
static
void
add_exemplar_alignment(const alignment& al,
                       const unsigned max_indel_size,
                       const bool is_remove_leading_edge_indels,
                       const bool is_remove_trailing_edge_indels,
                       std::vector<alignment>& exal)
{
    if (! al.is_realignable(max_indel_size)) return;

#ifdef DEBUG_ALIGN
    log_os << "Adding Exemplar Alignment."
           << " remove leading: " << is_remove_leading_edge_indels
           << " remove trailing: " << is_remove_trailing_edge_indels
           << "\n";
#endif

    // right now we:
    // (1) optionally remove soft-clipping and force the clipped regions to match
    // (2) remove edge indels for grouper reads
    // TODO -- find a better way to handle this case:
    //
    const alignment* al_ptr(&al);
    alignment nial;
    if ((is_remove_leading_edge_indels || is_remove_trailing_edge_indels) & is_edge_readref_len_segment(al.path))
    {
        nial=matchify_edge_indels(*al_ptr,is_remove_leading_edge_indels,is_remove_trailing_edge_indels);
        al_ptr=&nial;
    }

    alignment nscal;
    if (is_soft_clipped(al.path))
    {
        nscal=matchify_edge_soft_clip(*al_ptr);
        al_ptr=&nscal;
    }

    // check that this candidate exemplar does not already exist:
    for (alignment& exemplar : exal)
    {
        if (check_and_adjust_exemplar(*al_ptr,exemplar)) return;
    }

    // no compatible alignment found! Add this alignment as a new exemplar:
    exal.push_back(*al_ptr);

#ifdef DEBUG_ALIGN
    log_os << "Added exemplar " << *al_ptr << "\n";
#endif

}



// reduce discovery alignments to a set of non-redundant exemplars:
//
void
get_exemplar_alignments(const starling_base_options& opt,
                        const read_segment& rseg,
                        std::vector<alignment>& exal)
{
    exal.clear();

    // get exemplar from read mapper:
    const alignment& al(rseg.genome_align());
    if (! al.empty())
    {
        const std::pair<bool,bool> end_pin(rseg.get_segment_edge_pin());
        const bool is_remove_leading_edge_indels(! end_pin.first);
        const bool is_remove_trailing_edge_indels(! end_pin.second);
        add_exemplar_alignment(al,opt.max_indel_size,
                               is_remove_leading_edge_indels,
                               is_remove_trailing_edge_indels,
                               exal);
    }

#ifdef DEBUG_ALIGN
    log_os << "VARMIT: Final exemplar set:\n";
    for (const alignment& exemplar : exal)
    {
        log_os << "exemplar: " << exemplar;
    }
#endif
}



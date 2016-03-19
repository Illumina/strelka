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

#include "SequenceErrorCountsPosProcessor.hh"
#include "blt_common/ref_context.hh"
#include "common/OutStream.hh"
#include "starling_common/starling_indel_call_pprob_digt.hh"
#include "starling_common/starling_diploid_indel.hh"
#include "starling_common/starling_indel_error_prob.hh"



SequenceErrorCountsPosProcessor::
SequenceErrorCountsPosProcessor(
    const SequenceErrorCountsOptions& opt,
    const SequenceErrorCountsDerivOptions& dopt,
    const reference_contig_segment& ref,
    const SequenceErrorCountsStreams& streams)
    : base_t(opt,dopt,ref,streams,1),
      _opt(opt),
      _dopt(dopt),
      _streams(streams)
{
    // check that we have write permission on the output file early:
    {
        OutStream outs(opt.countsFilename);
    }

    // setup indel syncronizers:
    {
        sample_info& normal_sif(sample(0));

        if (dopt.is_max_depth())
        {
            if (opt.max_candidate_indel_depth_factor > 0.)
            {
                _max_candidate_normal_sample_depth = (opt.max_candidate_indel_depth_factor * dopt.max_depth);
            }
        }

        if (opt.max_candidate_indel_depth > 0.)
        {
            if (_max_candidate_normal_sample_depth > 0.)
            {
                _max_candidate_normal_sample_depth = std::min(_max_candidate_normal_sample_depth,static_cast<double>(opt.max_candidate_indel_depth));
            }
            else
            {
                _max_candidate_normal_sample_depth = opt.max_candidate_indel_depth;
            }
        }

        indel_sync_data isdata;
        isdata.register_sample(normal_sif.indel_buff,normal_sif.estdepth_buff,normal_sif.estdepth_buff_tier2,
                               normal_sif.sample_opt, _max_candidate_normal_sample_depth, 0);
        normal_sif.indel_sync_ptr.reset(new indel_synchronizer(opt, ref, dopt.countCache, isdata, 0));
    }
}

SequenceErrorCountsPosProcessor::
~SequenceErrorCountsPosProcessor()
{
    _counts.save(_opt.countsFilename.c_str());
}



void
SequenceErrorCountsPosProcessor::
reset()
{
    base_t::reset();
}



void
SequenceErrorCountsPosProcessor::
process_pos_indel_single_sample_digt(
    const pos_t pos,
    const unsigned sample_no)
{
    static const unsigned maxHpolLength(20);

    // note multi-sample status -- can still be called only for one sample
    // and only for sample 0. working on generalization:
    //
    if (sample_no!=0) return;



    // Current multiploid indel model can handle a het or hom indel
    // allele vs. reference, or two intersecting non-reference indel
    // alleles. (note that indel intersection is evaluated only in
    // terms of breakpoints -- so, for instance, a small het deletion
    // could occur within a large het deletion and the two would be
    // treated as non-interacting -- this is just an artifact of how
    // the methods are coded,)
    //

    typedef indel_buffer::const_iterator ciiter;

    sample_info& sif(sample(sample_no));


    // candidate indles are turned off above a certain depth relative to chromosome
    // mean, so we make sure to turn off counts as well to prevent a reference bias:
    {
        if (_max_candidate_normal_sample_depth > 0.)
        {
            const unsigned estdepth(sif.estdepth_buff.val(pos-1));
            const unsigned estdepth2(sif.estdepth_buff_tier2.val(pos-1));
            if ((estdepth+estdepth2) > _max_candidate_normal_sample_depth) return;
        }
    }

    ciiter it(sif.indel_buff.pos_iter(pos));
    const ciiter it_end(sif.indel_buff.pos_iter(pos+1));

    std::set<SequenceErrorContext> representedContexts;

    for (; it!=it_end; ++it)
    {
        const indel_key& ik(it->first);
        const indel_data& id(get_indel_data(it));
        const bool isForcedOutput(id.is_forced_output);

        if (! (ik.type == INDEL::DELETE || ik.type == INDEL::INSERT)) continue;

        if (! isForcedOutput)
        {
            if (! sif.indel_sync().is_candidate_indel(ik,id)) continue;
        }

        // TODO implement indel overlap resolution
        //
        // punt conflict resolution for now....
        {
            // indel_report_info needs to be run first now so that
            // local small repeat info is available to the indel
            // caller

            // sample-independent info:
            starling_indel_report_info iri;
            get_starling_indel_report_info(ik,id,_ref,iri);

            // STARKA-248 filter invalid indel
            /// TODO: filter this issue earlier (occurs as, e.g. 1D1I which matches ref)
            if (iri.vcf_indel_seq == iri.vcf_ref_seq) continue;

            double indel_error_prob(0);
            double ref_error_prob(0);

            static const bool is_tier2_pass(false);
            static const bool is_use_alt_indel(true);

            starling_diploid_indel dindel;
            dindel.is_forced_output = isForcedOutput;
            dindel.is_zero_coverage = false;

            {
                // check whether we're in a haploid/noploid region, for indels just check
                // start position and end position, approximating that the whole
                // region in between has the same ploidy, for any anomalous state
                // revert to 'noploid':
                const int indelLeftPloidy(get_ploidy(ik.pos));
                const int indelRightPloidy(get_ploidy(ik.right_pos()));

                if (indelLeftPloidy == indelRightPloidy)
                {
                    dindel.ploidy = indelLeftPloidy;
                }
                else
                {
                    dindel.ploidy = 0;
                }
            }

            _dopt.incaller().starling_indel_call_pprob_digt(
                _opt,_dopt,
                sif.sample_opt,
                indel_error_prob,ref_error_prob,
                ik,id,is_use_alt_indel,dindel);

            // sample-specific info: (division doesn't really matter
            // in single-sample case)
            starling_indel_sample_report_info isri;
            get_starling_indel_sample_report_info(_dopt,ik,id,sif.bc_buff,
                                                  is_tier2_pass,is_use_alt_indel,isri);

            SequenceErrorContext context;
            if (ik.type == INDEL::DELETE)
            {
                context.itype = INDEL_TYPE::DELETE;
            }
            else if(ik.type == INDEL::INSERT)
            {
                context.itype = INDEL_TYPE::INSERT;
            }
            else
            {
                assert(false);
            }

            context.hpolLength = 1;
            if ((iri.repeat_unit_length==1) && (iri.ref_repeat_count>1))
            {
                context.hpolLength = std::min(maxHpolLength, iri.ref_repeat_count);
            }

            // this would occur for, say, two different types of deletion noise at a homopolymer
            // that we don't currently pull apart -- really we should be summing all of these
            // types together - this is an approximation to get by....
            if (representedContexts.count(context)) continue;

            representedContexts.insert(context);

            SequenceErrorContextObservation obs;
            obs.signalCount = isri.n_q30_indel_reads;
            obs.refCount = (isri.n_q30_ref_reads + isri.n_q30_alt_reads);

            _counts.addError(context,obs,isri.depth);
        }
    }

    // add all the backgrounds that haven't been covered already
    {
        // background depth is always one minus position to be consistent with indel report:
        const pos_t depth_pos(pos-1);
        const snp_pos_info& spi(sif.bc_buff.get_pos(depth_pos));
        const unsigned depth(spi.calls.size());

        unsigned leftHpolSize(get_left_shifted_hpol_size(pos,_ref));

        SequenceErrorContext context;
        SequenceBackgroundObservation obs;
        obs.depth = depth;

        for (unsigned indelType(0);indelType<INDEL_TYPE::SIZE;++indelType)
        {
            context.itype = static_cast<INDEL_TYPE::index_t>(indelType);

            // always add hpol=1:
            context.hpolLength = 1;
            if (! representedContexts.count(context))
            {
                _counts.addBackground(context,obs);
            }

            // also check for a more specific reference context
            if (leftHpolSize<=1) continue;
            context.hpolLength = std::min(maxHpolLength,leftHpolSize);
            if (! representedContexts.count(context))
            {
                _counts.addBackground(context,obs);
            }
        }
    }
}

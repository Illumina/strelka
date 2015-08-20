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
 *  Created on: Jun 4, 2015
 *      Author: jduddy
 */

#include "variant_prefilter_stage.hh"
#include "calibration_models.hh"

variant_prefilter_stage::variant_prefilter_stage(const calibration_models& model, variant_pipe_stage_base& destination)
    : variant_pipe_stage_base(destination)
    , _model(model)
{
}

static
void
set_site_gt(const diploid_genotype::result_set& rs,
            site_modifiers& smod)
{
    smod.max_gt=rs.max_gt;
    smod.gqx=rs.max_gt_qphred;
    smod.gq  = 2;
}


void variant_prefilter_stage::add_site_modifiers(
    const site_info& si,
    site_modifiers& smod,
    const calibration_models& model)
{
    smod.clear();
    smod.is_unknown=(si.ref=='N');
    smod.is_used_covered=(si.n_used_calls!=0);
    smod.is_covered=(si.smod.is_used_covered || si.n_unused_calls!=0);

    if     (smod.is_unknown)
    {
        smod.gqx=0;
        smod.gq=0;
        smod.max_gt=0;
    }
    else if (si.dgt.genome.max_gt != si.dgt.poly.max_gt)
    {
        smod.gqx=0;
        smod.gq=si.dgt.poly.max_gt_qphred;
        smod.max_gt=si.dgt.poly.max_gt;
    }
    else
    {
        if (si.dgt.genome.max_gt_qphred<si.dgt.poly.max_gt_qphred)
        {
            set_site_gt(si.dgt.genome,smod);
        }
        else
        {
            set_site_gt(si.dgt.poly,smod);
        }
        smod.gq=si.dgt.poly.max_gt_qphred;
    }

    model.clasify_site(si, smod);
}




void variant_prefilter_stage::process(site_info& si)
{
    add_site_modifiers(si, si.smod, _model);
    if (si.dgt.is_haploid())
    {
        if (si.smod.max_gt == si.dgt.ref_gt)
        {
            si.smod.modified_gt=MODIFIED_SITE_GT::ZERO;
        }
        else
        {
            si.smod.modified_gt=MODIFIED_SITE_GT::ONE;
        }
    }
    else if (si.dgt.is_noploid())
    {
        if (! si.is_print_unknowngt())
        {
            si.smod.set_filter(VCF_FILTERS::PloidyConflict);
        }
    }

    _sink->process(si);

}
void variant_prefilter_stage::process(indel_info& ii)
{
    // add filter for all indels in no-ploid regions:
    if (ii.dindel.is_noploid())
    {
        ii.imod().set_filter(VCF_FILTERS::PloidyConflict);
    }

    _sink->process(ii);

}


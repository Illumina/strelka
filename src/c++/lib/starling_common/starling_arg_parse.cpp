//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2017 Illumina, Inc.
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

#include "blt_common/blt_arg_parse_util.hh"
#include "starling_common/starling_arg_parse.hh"

//#define DEBUG_PARSER

#ifdef DEBUG_PARSER
#include "blt_util/log.hh"
#endif


void
legacy_starling_arg_parse(
    arg_data& ad,
    starling_base_options& opt)
{
    const prog_info& pinfo(ad.pinfo);

    if ((ad.size()==1) ||
        ((ad.size()==2) && ((ad.argstr[1] == "-h") || (ad.argstr[1] == "-help") ||
                            (ad.argstr[1] == "--help") || (ad.argstr[1] == "-")))) pinfo.usage();

    opt.cmdline = ad.cmdline;

    bool is_min_qscore_set(false);
    bool is_min_pascore_set(false);

    bool is_bsnp_ssd_no_mismatch(false);
    bool is_bsnp_ssd_one_mismatch(false);

    bool is_max_indel_size(false);

    bool is_inmp(false);

    bool local_is_min_vexp(false);

    const unsigned as(ad.size());
    for (unsigned i(0); i<as; ++i)
    {
        if (ad.argmark[i]) continue;

        if (ad.argstr[i]=="-bsnp-ssd-no-mismatch")
        {
            set_xrange_arg(i,ad,is_bsnp_ssd_no_mismatch,opt.bsnp_ssd_no_mismatch,true);
        }
        else if (ad.argstr[i]=="-bsnp-ssd-one-mismatch")
        {
            set_xrange_arg(i,ad,is_bsnp_ssd_one_mismatch,opt.bsnp_ssd_one_mismatch,true);
        }
        else if (ad.argstr[i]=="-bsnp-diploid-het-bias")
        {
            set_xrange_arg(i,ad,opt.is_bsnp_diploid_het_bias,opt.bsnp_diploid_het_bias,true);
        }
        else if (ad.argstr[i]=="-min-qscore")
        {
            set_arg(i,ad,is_min_qscore_set,opt.min_qscore);
        }
        else if (ad.argstr[i]=="-max-window-mismatch")
        {
            int max_win_mismatch_tmp;
            set_win_arg(i,ad,opt.is_max_win_mismatch,max_win_mismatch_tmp,opt.max_win_mismatch_flank_size);
            if (max_win_mismatch_tmp<0)
            {
                pinfo.usage("first argument following -max-window-mismatch must be a non-negative integer\n");
            }
            opt.max_win_mismatch=max_win_mismatch_tmp;
        }
        else if (ad.argstr[i]=="-min-mapping-quality")
        {
            set_arg(i,ad,is_min_pascore_set,opt.min_mapping_quality);
        }
        else if (ad.argstr[i]=="-indel-nonsite-match-prob")
        {
            set_xrange_arg(i,ad,is_inmp,opt.indel_nonsite_match_prob,true);
        }
        else if (ad.argstr[i]=="-min-vexp")
        {
            set_xrange_arg(i,ad,local_is_min_vexp,opt.min_vexp,true);
            opt.is_min_vexp=local_is_min_vexp;
        }
        else if (ad.argstr[i]=="-max-indel-size")
        {
            set_arg(i,ad,is_max_indel_size,opt.maxIndelSize);
        }
        else if (ad.argstr[i]=="-h")
        {
            pinfo.usage();
        }
        else
        {
            continue;
        }

        ad.argmark[i] = true;
    }

    ad.finalize_args();
}

//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2018 Illumina, Inc.
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


#include "somaticAlleleUtil.hh"
#include "somatic_call_shared.hh"
#include "blt_util/qscore.hh"



void
updateAlleleEVSScore(
    const VariantScoringModelServer& varModel,
    const result_set& rs,
    strelka_shared_modifiers& smod)
{
    smod.isEVS = true;
    smod.EVS = varModel.scoreVariant(smod.features.getAll());

    static const double maxEmpiricalVariantScore(60.0);
    smod.EVS = std::min(error_prob_to_phred(smod.EVS), maxEmpiricalVariantScore);

    if (rs.ntype != NTYPE::REF)
    {
        smod.EVS = 0.0;
    }
}

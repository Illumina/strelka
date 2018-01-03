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

///
/// \author Chris Saunders
///

#pragma once

#include "somaticVariantEmpiricalScoringFeatures.hh"
#include "calibration/VariantScoringModelServer.hh"

#include <cassert>

#include <array>
#include <bitset>
#include <iosfwd>


namespace SOMATIC_VARIANT_VCF_FILTERS
{

enum index_t
{
    // Shared:
    HighDepth,
    // SNV only:
    LowEVSsnv,
    BCNoise,
    SpanDel,
    QSS_ref,
    // indel only:
    LowEVSindel,
    IndelBCNoise,
    QSI_ref,
    LowDepth,
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
    case LowEVSsnv:
        return "LowEVS";
    case BCNoise:
        return "BCNoise";
    case SpanDel:
        return "SpanDel";
    case QSS_ref:
        return "QSS_ref";
    case IndelBCNoise:
        return "BCNoise";
    case QSI_ref:
        return "QSI_ref";
    case LowEVSindel:
        return "LowEVS";
    case LowDepth:
        return "LowDepth";
    default:
        assert(false && "Unknown vcf filter id");
        return nullptr;
    }
}
}


struct strelka_filter_keeper
{
    strelka_filter_keeper()
    {
        clear();
    }

    void
    set(const SOMATIC_VARIANT_VCF_FILTERS::index_t i)
    {
        if (_filters.test(i))
        {
            assert(false && "Setting filter twice");
        }
        _filters.set(i);
    }

    void
    write(
        std::ostream& os) const
    {
        if (_filters.none())
        {
            os << "PASS";
            return;
        }

        bool is_sep(false);
        for (unsigned i(0); i<SOMATIC_VARIANT_VCF_FILTERS::SIZE; ++i)
        {
            if (! _filters.test(i)) continue;

            if (is_sep)
            {
                os << ";";
            }
            else
            {
                is_sep=true;
            }
            os << SOMATIC_VARIANT_VCF_FILTERS::get_label(i);
        }
    }

    void
    clear()
    {
        _filters.reset();
    }

private:
    std::bitset<SOMATIC_VARIANT_VCF_FILTERS::SIZE> _filters;
};


struct strelka_shared_modifiers
{
    bool isEVS = false;
    double EVS = 0;
    strelka_filter_keeper filters;
    VariantScoringFeatureKeeper features;
    VariantScoringFeatureKeeper dfeatures;
protected:
    /// protected so that subclass is resposible for lifetime
    /// of featureSets
    strelka_shared_modifiers(
        const FeatureSet& featureSet,
        const FeatureSet& developmentFeatureSet)
        : features(featureSet),
          dfeatures(developmentFeatureSet)
    {}
};


struct strelka_shared_modifiers_snv : public strelka_shared_modifiers
{
    strelka_shared_modifiers_snv()
        : strelka_shared_modifiers(
              SOMATIC_SNV_SCORING_FEATURES::getInstance(),
              SOMATIC_SNV_SCORING_DEVELOPMENT_FEATURES::getInstance())
    {}
};


struct strelka_shared_modifiers_indel : public strelka_shared_modifiers
{
    strelka_shared_modifiers_indel()
        : strelka_shared_modifiers(
              SOMATIC_INDEL_SCORING_FEATURES::getInstance(),
              SOMATIC_INDEL_SCORING_DEVELOPMENT_FEATURES::getInstance())
    {}
};

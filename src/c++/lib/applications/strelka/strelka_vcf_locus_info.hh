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

///
/// \author Chris Saunders
///

#pragma once

#include "strelkaVQSRFeatures.hh"
#include "calibration/SerializedModel.hh"

#include <cassert>

#include <bitset>
#include <iosfwd>


namespace STRELKA_VCF_FILTERS
{

enum index_t
{
    // SNVs and indels:
    HighDepth,
    LowQscore,
    // SNVs only:
    BCNoise,
    SpanDel,
    QSS_ref,
    // indels only:
    Repeat,
    iHpol,
    IndelBCNoise,
    QSI_ref,
    Nonref,
    TOR,
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
        //return "DP"; // old strelka workflow name
    case BCNoise:
        return "BCNoise";
    case SpanDel:
        return "SpanDel";
    case QSS_ref:
        return "QSS_ref";
    case Repeat:
        return "Repeat";
    case iHpol:
        return "iHpol";
    case IndelBCNoise:
        return "BCNoise";
    case QSI_ref:
        return "QSI_ref";
    case LowQscore:
        return "LowQscore";
    case TOR:
        return "TOR";
    case Nonref:
        return "Nonref";
    default:
        assert(false && "Unknown vcf filter id");
        return nullptr;
    }
}
}

template<class _vqsr_featureset = STRELKA_SNV_VQSR_FEATURES>
struct strelka_shared_modifiers
{
    strelka_shared_modifiers()
    {
        clear();
    }

    void
    set_filter(const STRELKA_VCF_FILTERS::index_t i)
    {
        filters.set(i);
    }

    void
    set_feature(const typename _vqsr_featureset::index_t i,double val)
    {
        if (_isFeatureSet.test(i))
        {
            assert(false && "Set VQSR feature twice");
        }
        _featureVal[i] = val;
        _isFeatureSet.set(i);
    }

    double
    get_feature(const typename _vqsr_featureset::index_t i) const
    {
        assert(_isFeatureSet.test(i));
        return this->_featureVal.at(i);
    }

    bool
    test_feature(const typename _vqsr_featureset::index_t i) const
    {
        return _isFeatureSet[i];
    }

    const feature_type&
    get_features() const
    {
        return _featureVal;
    }

    inline
    void
    write_filters(
            std::ostream& os) const
    {
        if (filters.none())
        {
            os << "PASS";
            return;
        }

        bool is_sep(false);
        for (unsigned i(0); i<STRELKA_VCF_FILTERS::SIZE; ++i)
        {
            if (! filters.test(i)) continue;

            if (is_sep)
            {
                os << ";";
            }
            else
            {
                is_sep=true;
            }
            os << STRELKA_VCF_FILTERS::get_label(i);
        }
    }

    //For debug only
    void
    write_feature(
            std::ostream& os) const
    {
        os << "\n #FEAT ";
        for (auto it = _featureVal.cbegin(); it != _featureVal.cend(); ++it)
            os << _vqsr_featureset::get_feature_label(it->first) << "=" << it->second << "; ";
        os << "\n";
    }

    void
    clear()
    {
        filters.reset();
        _isFeatureSet.reset();
        _featureVal.clear();
        isQscore=false;
        Qscore = 0;
    }

    std::bitset<_vqsr_featureset::SIZE> filters;
    bool isQscore;
    double Qscore;

private:
    std::bitset<_vqsr_featureset::SIZE> _isFeatureSet;
    feature_type _featureVal; // holds VQSR features
};


template<class _vqsr_featureset>
std::ostream&
operator<<(
        std::ostream& os,
        const strelka_shared_modifiers<_vqsr_featureset>& shmod)
{
    os << " filters: ";
    shmod.write_filters(os);

    return os;
}

typedef strelka_shared_modifiers<STRELKA_SNV_VQSR_FEATURES> strelka_shared_modifiers_snv;
typedef strelka_shared_modifiers<STRELKA_INDEL_VQSR_FEATURES> strelka_shared_modifiers_indel;


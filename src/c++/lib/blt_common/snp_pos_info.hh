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

#include "blt_common/hapscore.hh"
#include "blt_common/MapqTracker.hh"
#include "blt_util/fastRanksum.hh"
#include "blt_util/MeanTracker.hh"
#include "blt_util/qscore.hh"
#include "blt_util/seq_util.hh"

#include <cstdint>

#include <iosfwd>
#include <vector>

//#define BC_DEBUG

///mothballing this to support VS2013
#if 0
/// max 0-indexed value available with bitCount bits
inline
constexpr
unsigned
bitMaxval(const unsigned bitCount)
{
    return ((1 << bitCount) - 1);
}
#endif


struct base_call
{
    base_call(const uint8_t init_base_id,
              const uint8_t init_qscore,
              const bool init_ifs,
#ifdef BC_DEBUG
              const uint16_t init_read_pos,
              const uint16_t init_read_size,
#else
              const uint16_t,
              const uint16_t,
#endif
              const bool init_is_call_filter,
              const bool init_is_neighbor_mismatch,
              const bool init_is_tscf)
        : qscore(std::min(static_cast<unsigned>(init_qscore),qscore_max)),
#ifdef BC_DEBUG
          read_pos(init_read_pos), read_size(init_read_size),
#endif
          base_id(init_base_id),
          is_fwd_strand(init_ifs),
          is_neighbor_mismatch(init_is_neighbor_mismatch),
          is_call_filter(init_is_call_filter),
          is_tier_specific_call_filter(init_is_tscf)
    {
        qphred_cache::qscore_check(qscore,"basecall quality");
    }

    // pull quality value transformations from caching functions:
    double
    error_prob() const
    {
        return qphred_to_error_prob(static_cast<int>(qscore));
    }

    double
    ln_error_prob() const
    {
        return qphred_to_ln_error_prob(static_cast<int>(qscore));
    }

    double
    ln_comp_error_prob() const
    {
        return qphred_to_ln_comp_error_prob(static_cast<int>(qscore));
    }

    uint16_t
    get_qscore() const
    {
        return qscore;
    }

private:
    static const unsigned qscore_max;

    enum { qscore_bits = 6 };

uint16_t qscore:
    qscore_bits;
public:
    uint16_t base_id:4;
    uint16_t is_fwd_strand:1;
    uint16_t is_neighbor_mismatch:1;
    uint16_t is_call_filter:1; // filtered from snp-calling
    uint16_t is_tier_specific_call_filter:1;
#ifdef BC_DEBUG
    uint16_t read_pos; // zero-indexed cycle number
    uint16_t read_size;
#endif
};

std::ostream& operator<<(std::ostream& os,const base_call& bc);


/// \brief Captures basecall information for a single position 'pileup'
struct snp_pos_info
{
    snp_pos_info()
    {
        clear();
    }

    void
    clear()
    {
        _is_ref_set=false;
        _ref_base='N';
        is_n_ref_warn=false;
        calls.clear();
        tier2_calls.clear();
        spanningDeletionReadCount=0;
        submappedReadCount=0;
        mapqTracker.clear();
        hap_set.clear();
        mq_ranksum.clear();
        baseq_ranksum.clear();
        readPositionRankSum.clear();
        distanceFromReadEdge.clear();
        altAlleleReadPositionInfo.clear();

        spanningIndelPloidyModification = 0;
    }

    /// \brief Get basecall counts from the pileup
    ///
    /// Basecalls are not counted if they are unknown or have quality score less than \p min_qscore
    template <typename T>
    void
    getBasecallCounts(
        std::array<T, N_BASE>& base_count,
        const int min_qscore = 0) const
    {
        for (unsigned i(0); i<N_BASE; ++i) base_count[i] = 0;

        for (const auto& call : calls)
        {
            if (call.base_id==BASE_ID::ANY) continue;
            if (call.get_qscore()<min_qscore) continue;
            base_count[call.base_id]++;
        }
    }

    unsigned
    get_most_frequent_alt_id(const unsigned ref_gt) const
    {
        unsigned alt_count[BASE_ID::SIZE] = {};
        for (const base_call& tbc : calls)
        {
            const uint8_t obs_id(tbc.base_id);
            if (obs_id == ref_gt || obs_id == BASE_ID::ANY) continue;
            ++alt_count[obs_id];
        }
        unsigned alt_id = ref_gt;
        unsigned max_count = 0;
        for (unsigned base_id(0); base_id<BASE_ID::SIZE; ++base_id)
        {
            if (alt_count[base_id] > max_count)
            {
                if (base_id == ref_gt) continue;
                max_count = alt_count[base_id];
                alt_id = base_id;
            }
        }

        return alt_id;
    }

    bool
    is_ref_set() const
    {
        return _is_ref_set;
    }

    void
    set_ref_base(char base)
    {
        _is_ref_set = true;
        _ref_base = base;
    }

    char
    get_ref_base() const
    {
        assert(_is_ref_set);
        return _ref_base;
    }

    /// \returns the read-position rank sum
    double
    get_read_pos_ranksum() const;

    /// \returns the mapQ rank sum
    double
    get_mq_ranksum() const;

    /// \return the baseq rank sum
    double
    get_baseq_ranksum() const;

    double
    get_raw_baseQ() const;

    /// \return the raw pos sum
    double
    get_raw_pos() const;

    void
    print_known_counts(std::ostream& os,
                       const int min_qscore = 0) const;

    void
    print_known_qscore(std::ostream& os,
                       const int min_qscore = 0) const;

private:
    bool _is_ref_set;
    char _ref_base; // always fwd-strand base
public:
    bool is_n_ref_warn;
    std::vector<base_call> calls;
    std::vector<base_call> tier2_calls; // call not passing stringent quality criteria

    /// number of spanning deletion reads crossing the site
    unsigned spanningDeletionReadCount;

    /// number of submapped reads crossing the site
    ///
    /// note this could be usable,filtered or spanning deletion,
    /// all submapped reads get counted here:
    unsigned submappedReadCount;

    MapqTracker mapqTracker;

    mutable hap_set_t hap_set;

    //for calculating various rank-sum statistics
    fastRanksum mq_ranksum;
    fastRanksum baseq_ranksum;
    fastRanksum readPositionRankSum;

    /// Track summary stats for the distance of the variant locus from the edge of the read.
    ///
    /// This is used to support an RNA-Seq EVS feature.
    MeanTracker distanceFromReadEdge;

    /// Utility to store read position and total read length.
    struct ReadPositionInfo
    {
        uint16_t readPos;
        uint16_t readLength;
    };

    /// Read position of all non-reference allele observations.
    ///
    /// This is used to compute an allele position bias features in the somatic model.
    std::vector<ReadPositionInfo> altAlleleReadPositionInfo;

    int spanningIndelPloidyModification = 0;
};

std::ostream& operator<<(std::ostream& os,const snp_pos_info& pci);



struct extended_pos_info
{
    extended_pos_info(const snp_pos_info& pi_init,
                      const std::vector<float>& de_init)
        : pi(pi_init)
        , de(de_init) {}

    const snp_pos_info& pi;
    const std::vector<float>& de;
};

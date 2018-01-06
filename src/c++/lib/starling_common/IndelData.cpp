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

#include "IndelData.hh"
#include "AlleleReportInfoUtil.hh"

#include "blt_util/log.hh"
#include "blt_util/prob_util.hh"
#include "calibration/IndelErrorModel.hh"
#include "common/Exceptions.hh"

#include <iostream>
#include <sstream>


//#define DEBUG_ID



void
ReadPathScores::
insertAlt(
    const IndelKey& indelKey,
    const score_t a)
{
    const unsigned ais(static_cast<unsigned>(alt_indel.size()));
    if (ais < 2)
    {
        alt_indel.push_back(std::make_pair(indelKey,a));
    }
    else
    {
        unsigned min_index(ais);
        score_t min(a);
        for (unsigned i(0); i<ais; ++i)
        {
            if (alt_indel[i].second < min)
            {
                min = alt_indel[i].second;
                min_index = i;
            }
        }
        if (min_index<ais)
        {
            alt_indel[min_index] = std::make_pair(indelKey,a);
        }
    }
}



void
IndelSampleData::
addIndelObservation(
    const IndelObservationData& obs_data)
{
    const bool isAbstractObservation(obs_data.is_external_candidate || obs_data.is_forced_output);

    if (! isAbstractObservation)
    {
        using namespace INDEL_ALIGN_TYPE;

        evidence_t* insertTarget(nullptr);

        if (obs_data.is_noise)
        {
            insertTarget=&noise_read_ids;
        }
        else if (obs_data.iat == GENOME_TIER1_READ)
        {
            insertTarget=&tier1_map_read_ids;
        }
        else if (obs_data.iat == GENOME_TIER2_READ)
        {
            insertTarget=&tier2_map_read_ids;
        }
        else if (obs_data.iat == GENOME_SUBMAP_READ)
        {
            insertTarget=&submap_read_ids;
        }
        else
        {
            assert(false && "Unknown indel alignment type");
        }

        assert (insertTarget != nullptr);
        insertTarget->insert(obs_data.id);

        // assertion is removed because the same read can be inserted twice in active regions
        // TODO: revisit this later
//        const auto retval = insertTarget->insert(obs_data.id);
//        assert (retval.second);
    }
}



static
void
scaleIndelErrorRate(
    const double logScaleFactor,
    double& indelErrorRate)
{
    static const double minIndelErrorProb(0.0);
    static const double maxIndelErrorProb(0.5);

    indelErrorRate = std::min(indelErrorRate, maxIndelErrorProb);
    indelErrorRate = softMaxInverseTransform(indelErrorRate, minIndelErrorProb, maxIndelErrorProb);
    indelErrorRate += logScaleFactor;
    indelErrorRate = softMaxTransform(indelErrorRate, minIndelErrorProb, maxIndelErrorProb);
}



void
IndelSampleData::
initializeAuxInfo(
    const starling_base_options& opt,
    const starling_base_deriv_options& dopt,
    const unsigned sampleIndex,
    const IndelKey& indelKey,
    const AlleleReportInfo& indelReportInfo)
{
    double refToIndelErrorProb;
    double indelToRefErrorProb;
    dopt.getIndelErrorModel().getIndelErrorRate(sampleIndex, indelKey, indelReportInfo, refToIndelErrorProb, indelToRefErrorProb);

    if (opt.isIndelRefErrorFactor)
    {
        scaleIndelErrorRate(dopt.logIndelRefErrorFactor, indelToRefErrorProb);
    }

    _errorRates.refToIndelErrorProb.updateValue(refToIndelErrorProb);
    _errorRates.indelToRefErrorProb.updateValue(indelToRefErrorProb);

    dopt.getIndelErrorModel().getIndelErrorRate(sampleIndex, indelKey, indelReportInfo, refToIndelErrorProb, indelToRefErrorProb, true);

    _errorRates.candidateRefToIndelErrorProb.updateValue(refToIndelErrorProb);
    _errorRates.candidateIndelToRefErrorProb.updateValue(indelToRefErrorProb);
}



void
IndelData::
initializeAuxInfo(
    const starling_base_options& opt,
    const starling_base_deriv_options& dopt,
    const reference_contig_segment& ref)
{
    getAlleleReportInfo(_indelKey, ref, _reportInfo);

    const unsigned sampleCount(getSampleCount());
    for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
    {
        getSampleData(sampleIndex).initializeAuxInfo(opt, dopt, sampleIndex, _indelKey, _reportInfo);
    }
}



void
IndelData::
addIndelObservation(
    const unsigned sampleIndex,
    const IndelObservationData& observationData)
{
#ifdef DEBUG_ID
    log_os << "KATTER: input obs: " << observationData;
#endif

    // never reset the flags to false if they are true already
    if (! is_external_candidate) is_external_candidate=observationData.is_external_candidate;
    if (! isForcedOutput) isForcedOutput=observationData.is_forced_output;

    if (_indelKey.is_breakpoint())
    {
        _breakpointInsertSeq.addObservation(observationData.breakpointInsertionSequence);
    }
    else
    {
        if (not observationData.breakpointInsertionSequence.empty())
        {
            using namespace illumina::common;

            std::ostringstream oss;
            oss << "Indel observation sets breakpoint insertion sequence for a non-breakpoint allele: " << _indelKey << "\n";
            oss << "\tobservationData: " << observationData << "\n";
            BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
        }

    }

    getSampleData(sampleIndex).addIndelObservation(observationData);
}



std::ostream&
operator<<(
    std::ostream& os,
    const IndelObservationData& obs)
{
    os << "is_noise: " << obs.is_noise << "\n";
    os << "is_external: " << obs.is_external_candidate << "\n";
    os << "is_forced_output: " << obs.is_forced_output << "\n";
    os << "type: " << INDEL_ALIGN_TYPE::label(obs.iat) << "\n";
    os << "align_id: " << obs.id << "\n";
    os << "breakpointInsertionSequence: " << obs.breakpointInsertionSequence << "\n";
    return os;
}



static
void
report_indel_evidence_set(
    const IndelSampleData::evidence_t& e,
    const char* label,
    std::ostream& os)
{
    typedef IndelSampleData::evidence_t::const_iterator viter;
    viter i(e.begin()),i_end(e.end());
    for (unsigned n(0); i!=i_end; ++i)
    {
        os << label << " no: " << ++n << " id: " << *i << "\n";
    }
}





std::ostream&
operator<<(
    std::ostream& os,
    const ReadPathScores& rps)
{
    os << "ref: " << rps.ref
       << " indel: " << rps.indel
       << " nonAmbiguousBasesInRead: " << rps.nonAmbiguousBasesInRead;

    for (const auto& indel : rps.alt_indel)
    {
        const IndelKey& indelKey(indel.first);
        os << " alt-score: " << indel.second << " key: " << indelKey;
    }

    os << " ist1?: " << rps.is_tier1_read;

    return os;
}



void
BreakpointInsertSequenceManager::
_exception(const char* msg) const
{
    using namespace illumina::common;
    std::ostringstream oss;
    oss << "Exception in BreakpointInsertSequenceManager: " << msg;
    BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
}



void
BreakpointInsertSequenceManager::
_finalize()
{
    unsigned count(0);
    std::string& candidate(_consensus_seq);

    for (const auto& val : _obs)
    {
        if (val.first.size() < candidate.size()) continue;
        if (val.first.size() == candidate.size())
        {
            if (val.second <= count) continue;
        }
        candidate = val.first;
        count = val.second;
    }
    _consensus_seq = candidate;
    _obs.clear();
    _is_consensus=true;
}



std::ostream&
operator<<(
    std::ostream& os,
    const IndelSampleData& indelSampleData)
{
    report_indel_evidence_set(indelSampleData.tier1_map_read_ids,"tier1_map_read",os);
    report_indel_evidence_set(indelSampleData.tier2_map_read_ids,"tier2_map_read",os);
    report_indel_evidence_set(indelSampleData.submap_read_ids,"submap_read",os);
    report_indel_evidence_set(indelSampleData.noise_read_ids,"noise_read",os);

    {
        unsigned n(0);
        for (const auto& readLnp : indelSampleData.read_path_lnp)
        {
            os << "read_path_lnp no: " << ++n
               << " id: " << readLnp.first
               << " " << readLnp.second
               << "\n";
        }
    }

    report_indel_evidence_set(indelSampleData.suboverlap_tier1_read_ids,"suboverlap_tier1_read",os);
    report_indel_evidence_set(indelSampleData.suboverlap_tier2_read_ids,"suboverlap_tier2_read",os);

    return os;
}



std::ostream&
operator<<(
    std::ostream& os,
    const IndelData& indelData)
{
    os << "IndelKey: " << indelData._indelKey << "\n";
    os << "is_external_candidate: " << indelData.is_external_candidate << "\n";
    os << "is_forced_output: " << indelData.isForcedOutput << "\n";

    os << "breakpointInsertionSequence: " << indelData.getBreakpointInsertSeq() << "\n";

    const unsigned sampleCount(indelData.getSampleCount());
    for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
    {
        os << "BEGIN sample: " << sampleIndex << "\n";
        os << indelData.getSampleData(sampleIndex);
        os << "BEGIN sample: " << sampleIndex << "\n";
    }
    return os;
}

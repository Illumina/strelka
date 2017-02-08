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
/*
 *
 *  Created on: Sep 10, 2013
 *  Author: Morten Kallberg
 */

#include "codon_phaser.hh"

#include <array>
#include <functional>
#include <vector>

//#define DEBUG_CODON



#ifdef DEBUG_CODON
#include "blt_util/log.hh"
#endif


void
Codon_phaser::
process(std::unique_ptr<GermlineSiteLocusInfo> locusPtr)
{
    if (not _opt.isUseCodonPhaser)
    {
        // bypass phasing
        _sink->process(std::move(locusPtr));
        return;
    }

    auto curActiveRegionId(locusPtr->getActiveRegionId());

    if (curActiveRegionId < 0)
    {
        // locus is outside of active region
        outputBuffer();
        _sink->process(std::move(locusPtr));
    }
    else if (curActiveRegionId != _activeRegionId)
    {
        // new active region
        outputBuffer();

        _activeRegionId = curActiveRegionId;
        addSiteLocusToBuffer(std::move(locusPtr));
    }
    else
    {
        // locus is within the existing active region
        addSiteLocusToBuffer(std::move(locusPtr));
    }
}

void
Codon_phaser::
process(std::unique_ptr<GermlineIndelLocusInfo> locusPtr)
{
    if (not _opt.isUseCodonPhaser)
    {
        // bypass phasing
        _sink->process(std::move(locusPtr));
        return;
    }

    auto curActiveRegionId(locusPtr->getActiveRegionId());

    if (curActiveRegionId < 0)
    {
        // locus is outside of active region
        outputBuffer();
        _sink->process(std::move(locusPtr));
    }
    else if (curActiveRegionId != _activeRegionId)
    {
        // new active region starts
        outputBuffer();
        _activeRegionId = curActiveRegionId;
        addIndelLocusToBuffer(std::move(locusPtr));
    }
    else
    {
        // locus is within the existing active region
        addIndelLocusToBuffer(std::move(locusPtr));
    }
}


void
Codon_phaser::
create_phased_record(
        const unsigned /*sampleIndex*/)
{

}


void
Codon_phaser::
outputBuffer()
{
    for (unsigned sampleId(0); sampleId < _sampleCount; ++sampleId)
        outputBuffer(sampleId);
}

void
Codon_phaser::
outputBuffer(unsigned sampleId)
{
    if (not isBuffer(sampleId)) return;

    auto& indelLocusBuffer(_sampleIndelLocusBuffer[sampleId]);
    auto& siteLocusBuffer(_sampleSiteLocusBuffer[sampleId]);

    if (indelLocusBuffer.empty())
    {
        for (auto& siteLocus : siteLocusBuffer)
        {
            _sink->process(std::move(siteLocus));
        }
    }
    else if (siteLocusBuffer.empty())
    {
        for (auto& indelLocus : indelLocusBuffer)
        {
            _sink->process(std::move(indelLocus));
        }
    }
    else
    {
        // mix of SNVs and indels
        unsigned siteBufferIndex(0u);
        unsigned indelBufferIndex(0u);
        while (siteBufferIndex < siteLocusBuffer.size() and indelBufferIndex < indelLocusBuffer.size())
        {
            // process SNVs and indels by the order of pos
            pos_t sitePos(siteLocusBuffer[siteBufferIndex]->pos);
            pos_t indelPos(indelLocusBuffer[indelBufferIndex]->pos);
            if (sitePos < indelPos)
                _sink->process(std::move(siteLocusBuffer[siteBufferIndex++]));
            else
                _sink->process(std::move(indelLocusBuffer[indelBufferIndex++]));
        }

        while (siteBufferIndex < siteLocusBuffer.size())
            _sink->process(std::move(siteLocusBuffer[siteBufferIndex++]));

        while (indelBufferIndex < indelLocusBuffer.size())
            _sink->process(std::move(indelLocusBuffer[indelBufferIndex++]));
    }

    siteLocusBuffer.clear();
    indelLocusBuffer.clear();
}



void
Codon_phaser::
flush_impl()
{
    if (_opt.isUseCodonPhaser)
    {
        outputBuffer();
    }
}

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
#pragma once
#include "gvcf_locus_info.hh"


class variant_pipe_stage
{
public:
    virtual void process(site_info& si) { if (_sink) _sink->process(si); }
    virtual void process(indel_info& ii) { if (_sink) _sink->process(ii); }
    virtual void flush()
    {
        if (_sink)
            _sink->flush();
    }

    variant_pipe_stage(variant_pipe_stage& sink) : _sink(&sink) {}
    variant_pipe_stage() : _sink(nullptr) {}


protected:



    variant_pipe_stage* _sink;
};


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


class variant_pipe_stage_base
{
public:
    virtual void process(std::unique_ptr<site_info> si) { if (_sink) _sink->process(std::move(si)); }
    virtual void process(std::unique_ptr<indel_info> ii) { if (_sink) _sink->process(std::move(ii)); }
    virtual void flush()
    {
        if (_sink)
            _sink->flush();
    }

    explicit variant_pipe_stage_base(std::shared_ptr<variant_pipe_stage_base> sink) : _sink(sink) {}

    virtual ~variant_pipe_stage_base() {}

protected:
    variant_pipe_stage_base() : _sink(nullptr) {}

    template <class TDerived, class TBase>
    static std::unique_ptr<TDerived> downcast(std::unique_ptr<TBase> basePtr)
    {
        if (typeid(*basePtr) == typeid(TDerived))
        {
            return std::unique_ptr<TDerived>(dynamic_cast<TDerived*>(basePtr.release()));
        }
        throw std::bad_cast();
    }




    std::shared_ptr<variant_pipe_stage_base> _sink;
};


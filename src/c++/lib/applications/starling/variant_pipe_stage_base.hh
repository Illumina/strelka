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
#pragma once
#include "gvcf_locus_info.hh"


class variant_pipe_stage_base
{
public:
    virtual void process(std::unique_ptr<site_info> si)
    {
        if (_sink) _sink->process(std::move(si));
    }
    virtual void process(std::unique_ptr<indel_info> ii)
    {
        if (_sink) _sink->process(std::move(ii));
    }
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


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

#include <cmath>


struct LogValuePair
{
    double
    getValue() const
    {
        return _value;
    }

    double
    getLogValue() const
    {
        return _logValue;
    }

    void
    updateValue(const double value)
    {
        _value = value;
        _logValue = std::log(_value);
    }

    void
    updateLogValue(const double logValue)
    {
        _logValue = logValue;
        _value = std::exp(_logValue);
    }

private:
    static const double _log0;
    double _value = 0.;
    double _logValue = _log0;
};

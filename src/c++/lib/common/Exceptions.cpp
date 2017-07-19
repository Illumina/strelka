//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2017 Illumina, Inc.
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

/**
 ** \file
 ** \brief Implementation of the common exception mechanism.
 **
 ** \author Come Raczy
 **/

#include "Exceptions.hh"

#include <cstring>
#include <cerrno>
#include <ctime>


namespace illumina
{
namespace common
{

ExceptionData::ExceptionData(int errorNumber, const std::string& message)
    : errorNumber_(errorNumber), message_(message)
{
}

std::string ExceptionData::getContext() const
{
    static const unsigned bufferSize(256);
    char timeBuffer[bufferSize];

    const time_t result(time(nullptr));
    strftime(timeBuffer, bufferSize, "%c", localtime(&result));

    return std::string(timeBuffer) + ": " + std::string(strerror(errorNumber_)) + ": " + boost::diagnostic_information(*this);
}

IoException::IoException(int errorNumber, const std::string& message)
    : std::ios_base::failure(message)
    , ExceptionData(errorNumber, message)
{
}

UnsupportedVersionException::UnsupportedVersionException(const std::string& message)
    : std::logic_error(message)
    , ExceptionData(EINVAL, message)
{
}

InvalidParameterException::InvalidParameterException(const std::string& message)
    : std::logic_error(message)
    , ExceptionData(EINVAL, message)
{
}

InvalidOptionException::InvalidOptionException(const std::string& message)
    : std::logic_error(message)
    , ExceptionData(EINVAL, message)
{
}

PreConditionException::PreConditionException(const std::string& message)
    : std::logic_error(message)
    , ExceptionData(EINVAL, message)
{
}

PostConditionException::PostConditionException(const std::string& message)
    : std::logic_error(message)
    , ExceptionData(EINVAL, message)
{
}

}
}

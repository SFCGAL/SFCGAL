// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/Exception.h"

#include <format>
#include <utility>

// NOLINTBEGIN(bugprone-throw-keyword-missing)
namespace SFCGAL {

Exception::Exception() noexcept : _message("unknown exception") {}

Exception::Exception(std::string          message,
                     std::source_location location) noexcept
    : _message(std::move(message))
{
  _diagnostic = std::format("{}:{} in {}", location.file_name(),
                            location.line(), location.function_name());
}

Exception::~Exception() noexcept = default;

auto
Exception::what() const noexcept -> const char *
{
  return _message.c_str();
}

auto
Exception::diagnostic() const noexcept -> std::string
{
  return _diagnostic;
}

// Definitions of constructors and destructors for derived classes

GeometryInvalidityException::GeometryInvalidityException(
    std::string const &message, std::source_location location) noexcept
    : Exception(message, location)
{
}

GeometryInvalidityException::~GeometryInvalidityException() noexcept = default;

NotImplementedException::NotImplementedException(
    std::string const &message, std::source_location location) noexcept
    : Exception(message, location)
{
}

NotImplementedException::~NotImplementedException() noexcept = default;

InappropriateGeometryException::InappropriateGeometryException(
    std::string const &message, std::source_location location) noexcept
    : Exception(message, location)
{
}

InappropriateGeometryException::~InappropriateGeometryException() noexcept =
    default;

NonFiniteValueException::NonFiniteValueException(
    std::string const &message, std::source_location location) noexcept
    : Exception(message, location)
{
}

NonFiniteValueException::~NonFiniteValueException() noexcept = default;

WktParseException::WktParseException(std::string const   &message,
                                     std::source_location location) noexcept
    : Exception(message, location)
{
}

WktParseException::~WktParseException() noexcept = default;

WkbParseException::WkbParseException(std::string const   &message,
                                     std::source_location location) noexcept
    : Exception(message, location)
{
}

WkbParseException::~WkbParseException() noexcept = default;

} // namespace SFCGAL
// NOLINTEND(bugprone-throw-keyword-missing)

// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/io/io_utils.h"
#include "SFCGAL/Exception.h"
#include <cctype>
#include <cmath>

namespace SFCGAL::io::detail {

void
skipWhitespace(std::istream &in)
{
  while (in && std::isspace(static_cast<unsigned char>(in.peek()))) {
    in.get();
  }
}

auto
parseDouble(std::istream &in, const std::string &context) -> double
{
  double value = 0.0;
  if (!(in >> value)) {
    BOOST_THROW_EXCEPTION(
        Exception("Failed to parse floating-point value in " + context));
  }

  if (std::isnan(value)) {
    BOOST_THROW_EXCEPTION(Exception("NaN value encountered in " + context));
  }
  if (std::isinf(value)) {
    BOOST_THROW_EXCEPTION(
        Exception("Infinite value encountered in " + context));
  }

  return value;
}

auto
caseInsensitiveEqual(const std::string &a, const std::string &b) -> bool
{
  if (a.size() != b.size()) {
    return false;
  }
  for (size_t i = 0; i < a.size(); ++i) {
    if (std::tolower(static_cast<unsigned char>(a[i])) !=
        std::tolower(static_cast<unsigned char>(b[i]))) {
      return false;
    }
  }
  return true;
}

} // namespace SFCGAL::io::detail

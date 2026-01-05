// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_IO_UTILS_H_
#define SFCGAL_IO_UTILS_H_

#include <istream>
#include <string>

namespace SFCGAL::io::detail {

/**
 * @brief Skip whitespace characters in stream
 * @param[in,out] in The input stream
 */
void
skipWhitespace(std::istream &in);

/**
 * @brief Parse a double value with robustness checks for NaN and Inf
 * @param[in] in The input stream
 * @param[in] context Description of what is being parsed (for error messages)
 * @return The parsed double value
 * @throws SFCGAL::Exception If parsing fails or value is NaN/Inf
 */
auto
parseDouble(std::istream &in, const std::string &context) -> double;

/**
 * @brief Case-insensitive string comparison
 * @param[in] a First string
 * @param[in] b Second string
 * @return true if strings are equal (case-insensitive)
 */
auto
caseInsensitiveEqual(const std::string &a, const std::string &b) -> bool;

} // namespace SFCGAL::io::detail

#endif // SFCGAL_IO_UTILS_H_

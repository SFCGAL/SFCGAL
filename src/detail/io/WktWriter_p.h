// Copyright (c) 2026-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_IO_WKBWRITER_P_H_
#define SFCGAL_IO_WKBWRITER_P_H_

#include <cmath>

// ----------------------------------------------------------------------------------
// -- private interface
// ----------------------------------------------------------------------------------
/// @cond PRIVATE_SECTION

namespace SFCGAL::detail::io {

static auto
fixZeroNeg(double val, int precision) -> double
{
  if (std::abs(val) < std::pow(10, -precision)) {
    return 0;
  }
  return val;
}

static auto
fixZeroNegForWeights(double val, int precision) -> double
{
  // NURBS weights must be positive
  if (val <= 0) {
    // Convert non-positive weights to minimum positive value
    return std::pow(10, -precision);
  }

  // For positive weights, preserve small values instead of converting to 0
  if (val < std::pow(10, -precision)) {
    return std::pow(10, -precision);
  }

  return val;
}

} // namespace SFCGAL::detail::io

/// @endcond

#endif // SFCGAL_IO_WKBWRITER_P_H_

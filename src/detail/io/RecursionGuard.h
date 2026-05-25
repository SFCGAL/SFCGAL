// Copyright (c) 2025-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_DETAIL_IO_RECURSIONGUARD_H_
#define SFCGAL_DETAIL_IO_RECURSIONGUARD_H_

namespace SFCGAL::detail::io {

/**
 * Guard to increment/decrement a recursion depth counter.
 *
 * Used in parsers and exporters (WKB, WKT, OBJ, STL, VTK) to guard
 * against unbounded recursion depth.  The depth counter is incremented
 * on construction and decremented on destruction, even if an exception
 * is thrown.
 */
struct RecursionGuard {
  /// Reference to the recursion depth counter being guarded
  int &depth;

  /**
   * @brief Increments the recursion depth counter
   * @param depthRef reference to the recursion depth counter
   */
  explicit RecursionGuard(int &depthRef) : depth(depthRef) { ++depth; }

  RecursionGuard(const RecursionGuard &) = delete;
  RecursionGuard(RecursionGuard &&)      = delete;
  auto
  operator=(const RecursionGuard &) -> RecursionGuard & = delete;
  auto
  operator=(RecursionGuard &&) -> RecursionGuard & = delete;

  ~RecursionGuard() { --depth; }
};

} // namespace SFCGAL::detail::io

#endif

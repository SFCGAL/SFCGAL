// Copyright (c) 2026-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGALOP_OPERATIONS_REGISTRY_HPP
#define SFCGALOP_OPERATIONS_REGISTRY_HPP

#pragma once

#include "operations.hpp"
#include <vector>

/**
 * Registry providing access to all available geometric operations.
 *
 * This class exposes a static interface to retrieve the full list of
 * registered operations used by sfcgalop.
 */
struct OperationRegistry {
  /**
   * Get the list of all registered operations.
   *
   * @return Constant reference to the vector containing all operations.
   */
  static auto
  all() -> const std::vector<Operation> &;
};

#endif // SFCGALOP_OPERATIONS_REGISTRY_HPP_

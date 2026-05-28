// Copyright (c) 2026-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGALOP_OPERATIONS_PROCESSING_HPP
#define SFCGALOP_OPERATIONS_PROCESSING_HPP

#include "operations.hpp"
#include <vector>

namespace Operations {

/**
 * @brief Returns the list of supported processing operations.
 * @return A reference to the static vector of processing operations.
 */
auto
get_operations_processing() -> const std::vector<Operation> &;

} // namespace Operations

#endif // SFCGALOP_OPERATIONS_PROCESSING_HPP

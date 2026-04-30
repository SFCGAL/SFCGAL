// Copyright (c) 2026-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "operations_registry.hpp"

#include "operations_analysis.hpp"
#include "operations_collections.hpp"
#include "operations_construction.hpp"
#include "operations_constructors.hpp"
#include "operations_conversions.hpp"
#include "operations_metrics.hpp"
#include "operations_predicates.hpp"
#include "operations_set.hpp"
#include "operations_transformations.hpp"

namespace {

auto
build_all_operations() -> std::vector<Operation>
{
  std::vector<Operation> ops;

  auto concat = [&](const std::vector<Operation> &operations) -> void {
    ops.insert(ops.end(), operations.begin(), operations.end());
  };

  concat(Operations::operations_analysis);
  concat(Operations::operations_collections);
  concat(Operations::operations_construction);
  concat(Operations::operations_constructors);
  concat(Operations::operations_conversions);
  concat(Operations::operations_metrics);
  concat(Operations::operations_predicates);
  concat(Operations::operations_set);
  concat(Operations::operations_transformations);

  return ops;
}

} // namespace

auto
OperationRegistry::all() -> const std::vector<Operation> &
{
  static const std::vector<Operation> ops = build_all_operations();
  return ops;
}

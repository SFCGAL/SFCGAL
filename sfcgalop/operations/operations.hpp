// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGALOP_OPERATIONS_HPP
#define SFCGALOP_OPERATIONS_HPP

#include "../io.hpp"
#include <SFCGAL/Geometry.h>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

namespace Operations {

// Type alias for operation metadata: (name, category, description, help,
// params)
using OperationInfo =
    std::tuple<std::string, std::string, std::string, std::string, std::string>;

auto
execute_operation(const std::string &op_name, const std::string &op_arg,
                  const SFCGAL::Geometry *geom_a,
                  const SFCGAL::Geometry *geom_b)
    -> std::optional<OperationResult>;

void
print_available_operations();
auto
print_operation_help(const char *name) -> bool;
auto
get_all_operations_info() -> std::vector<OperationInfo>;

auto
operation_requires_second_geometry(const std::string &operation_name) -> bool;

} // namespace Operations

// Provide non-namespaced aliases for backward compatibility and convenience.
// These allow calling execute_operation() without the Operations:: prefix.
// If namespace pollution is a concern, use the fully qualified names instead.
using Operations::execute_operation;
using Operations::print_available_operations;
using Operations::print_operation_help;

/**
 * @brief Structure defining a geometry operation for sfcgalop CLI
 *
 * Contains all metadata and implementation details for a single geometric
 * operation that can be executed via the command line interface.
 */
struct Operation {
  std::string name; ///< Operation name (e.g., "area", "intersection")
  std::string
      category; ///< Category for grouping (e.g., "Metrics", "Set Operations")
  std::string description; ///< Short description for operation list
  bool        requires_b;  ///< Whether operation requires a second geometry
  std::string param_help;  ///< Detailed help text with parameters and examples
  std::string input;       ///< Input specification (A, A,B, A,params)
  std::string output; ///< Output type (G=Geometry, D=Double, B=Boolean, T=Text)
  /// Function implementing the operation
  std::function<std::optional<OperationResult>(
      const std::string &, const SFCGAL::Geometry *, const SFCGAL::Geometry *)>
      func;
};

auto
parse_params(const std::string &str) -> std::map<std::string, double>;

auto
parse_boolean_param(const std::map<std::string, double> &params,
                    const std::string &key, const std::string &param_str,
                    bool default_val = false) -> bool;

auto
parse_double(const std::string &str, double default_val = 0.0) -> double;

auto
trim(const std::string &str) -> std::string;

#endif // SFCGALOP_OPERATIONS_HPP

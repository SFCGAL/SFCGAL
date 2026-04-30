// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "operations.hpp"
#include "operations_registry.hpp"

#include "../util.hpp"

#include <CGAL/number_utils.h>

#ifndef M_PI
  #define M_PI 3.14159265358979323846
#endif

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>

namespace {

using namespace SFCGAL::sfcgalop::util;

// Define the OperationMatchResult struct specifically for this implementation
struct OperationMatchResult {
  bool                     found;
  std::string              resolved_name;
  std::optional<Operation> operation;
};

auto
find_operation(std::string_view op_name) -> OperationMatchResult
{
  // Normalize the input operation name
  std::string            normalized_input = normalize_operation_name(op_name);
  std::vector<Operation> operations       = OperationRegistry::all();
  auto                   it               = std::find_if(
      operations.begin(), operations.end(), [&](const Operation &op) {
        return normalize_operation_name(op.name) == normalized_input;
      });

  if (it != operations.end()) {
    return {true, std::string(op_name), *it};
  }
  return {false, std::string(op_name), std::nullopt};
}

} // namespace

/**
 * @brief Parse a string to double, returning a fallback on failure.
 *
 * Attempts to convert the given string to a double using std::stod.
 * If conversion fails (invalid format, out-of-range, etc.), the provided
 * default_val is returned.
 *
 * @param str Input string to parse.
 * @param default_val Value to return if parsing fails (default: 0.0).
 * @return double Parsed double on success, otherwise default_val.
 */
auto
parse_double(const std::string &str, double default_val) -> double
{
  try {
    return std::stod(str);
  } catch (...) {
    return default_val;
  }
}

/**
 * @brief Parse a string as a boolean value.
 *
 * Accepts various boolean representations in a case-insensitive manner:
 * - "true", "t", "1" -> true
 * - "false", "f", "0" -> false
 *
 * @param str String to parse as boolean
 * @param default_val Default value to return if parsing fails
 * @return Boolean value or default_val if parsing fails
 */
auto
parse_boolean(const std::string &str, bool default_val) -> bool
{
  std::string lower_str = str;
  std::transform(lower_str.begin(), lower_str.end(), lower_str.begin(),
                 [](unsigned char chr) -> int { return std::tolower(chr); });

  if (lower_str == "true" || lower_str == "t" || lower_str == "1") {
    return true;
  }
  if (lower_str == "false" || lower_str == "f" || lower_str == "0") {
    return false;
  }
  return default_val;
}

/**
 * @brief Extract and parse a boolean parameter from the parameter map.
 *
 * This function handles both string-based boolean parameters and legacy
 * double-based parameters (where 0.0 = false, non-zero = true).
 *
 * @param params Parameter map from parse_params
 * @param key Parameter name to extract
 * @param param_str Original parameter string for string-based parsing
 * @param default_val Default value if parameter not found
 * @return Boolean value
 */

auto
parse_boolean_param(const std::map<std::string, double> &params,
                    const std::string &key, const std::string &param_str,
                    bool default_val) -> bool
{
  auto it = params.find(key);
  if (it == params.end()) {
    return default_val;
  }

  // Try to extract the original string value from param_str for proper boolean
  // parsing
  std::string search_key = key + "=";
  auto        key_pos    = param_str.find(search_key);
  if (key_pos != std::string::npos) {
    auto value_start = key_pos + search_key.length();
    auto value_end   = param_str.find(',', value_start);
    if (value_end == std::string::npos) {
      value_end = param_str.length();
    }
    std::string value_str =
        param_str.substr(value_start, value_end - value_start);
    value_str = trim(value_str);

    // If it looks like a string boolean, parse it as such
    if (value_str == "true" || value_str == "false" || value_str == "t" ||
        value_str == "f" || value_str == "T" || value_str == "F" ||
        value_str == "True" || value_str == "False" || value_str == "TRUE" ||
        value_str == "FALSE") {
      return parse_boolean(value_str, default_val);
    }
  }

  // Fallback to standard if it's not 0(.0) it's true
  return it->second != 0.0;
}

/**
 * @brief Parse a comma-separated list of key=value pairs into a map.
 *
 * Parses `str` for entries of the form `key=value` separated by commas and
 * returns a map from each key to its value parsed as a double. Entries that
 * do not contain an '=' character are ignored. Values are converted using
 * `parse_double`; any parse fallback behavior is handled by that function.
 *
 * @param str Input string containing comma-separated `key=value` pairs.
 * @return std::map<std::string,double> Mapping of keys to parsed double values.
 */
/**
 * Trim leading and trailing whitespace from a string
 */
auto
trim(const std::string &str) -> std::string
{
  auto start = str.begin();
  auto end   = str.end();

  // Find first non-whitespace
  while (start != end &&
         (std::isspace(static_cast<unsigned char>(*start)) != 0)) {
    ++start;
  }

  // Find last non-whitespace (from the end)
  while (start != end &&
         (std::isspace(static_cast<unsigned char>(*(end - 1))) != 0)) {
    --end;
  }

  return {start, end};
}

auto
parse_params(const std::string &str) -> std::map<std::string, double>
{
  std::map<std::string, double> params;
  std::istringstream            stream(str);
  std::string                   param;

  while (std::getline(stream, param, ',')) {
    auto eq_pos = param.find('=');
    if (eq_pos != std::string::npos) {
      auto name_raw  = param.substr(0, eq_pos);
      auto value_raw = param.substr(eq_pos + 1);

      // Trim whitespace from both key and value to prevent lookup failures
      auto name_trimmed  = trim(name_raw);
      auto value_trimmed = trim(value_raw);

      params[name_trimmed] = parse_double(value_trimmed);
    }
  }
  return params;
}

namespace Operations {

/**
 * @brief Execute a named geometry operation from the registry.
 *
 * Looks up the operation identified by op_name and, if found, invokes its
 * callable with op_arg and the optional geometry pointers geom_a and geom_b.
 * Exceptions thrown by the operation are caught and result in a nullopt return.
 *
 * @param op_name Name of the operation to execute (must match an entry in the
 * operations registry).
 * @param op_arg Optional argument string passed to the operation (parsed by the
 * operation as needed).
 * @param geom_a First geometry input (may be null if the operation does not
 * require it).
 * @param geom_b Second geometry input (may be null; many operations that
 * require a second geometry will return nullopt when this is absent).
 * @return std::optional<OperationResult> The operation's result if the
 * operation exists and completes successfully; std::nullopt if the operation
 * name is not found or if an exception occurs during execution.
 */
auto
execute_operation(const std::string &op_name, const std::string &op_arg,
                  const SFCGAL::Geometry *geom_a,
                  const SFCGAL::Geometry *geom_b)
    -> std::optional<OperationResult>
{
  auto match_result = find_operation(op_name);

  if (match_result.found) {
    // Check for null geometry A before calling operation, but allow constructor
    // operations
    bool is_constructor = match_result.operation->name.find("make_") == 0;
    if (geom_a == nullptr && !is_constructor) {
      return std::nullopt;
    }

    try {
      return match_result.operation->func(op_arg, geom_a, geom_b);
    } catch (const std::exception &e) {
      std::cerr << "Operation error: " << e.what() << "\n";
      return std::nullopt;
    }
  }

  return std::nullopt;
}

/**
 * @brief Prints a categorized list of all registered geometric operations.
 *
 * Outputs each operation name, short description, and a notice when an
 * operation requires two geometries. Operations are grouped by their category
 * and printed to stdout. Uses simple terminal color/formatting for readability.
 */
void
print_available_operations()
{
  std::map<std::string, std::vector<const Operation *>> by_category;

  for (const auto &operation : OperationRegistry::all()) {
    by_category[operation.category].push_back(&operation);
  }

  // Use simple formatting for now - can be enhanced with TextUI::Table later
  for (const auto &[category, operation_list] : by_category) {
    std::cout << "\n\033[1;36m" << category << ":\033[0m\n"; // Cyan bold
    for (const auto *operation : operation_list) {
      std::cout << "  \033[1m" << std::setw(25) << std::left << operation->name
                << "\033[0m"
                << " " << operation->description;
      if (operation->requires_b) {
        std::cout << " \033[33m(requires 2 geometries)\033[0m"; // Yellow
      }
      std::cout << "\n";
    }
  }
}

/**
 * @brief Print detailed help for a named operation.
 *
 * Looks up an operation by its name and writes its metadata (name, category,
 * description) to stdout. If the operation requires a second geometry this is
 * indicated. If the name is unknown an error message is written to stderr.
 *
 * @param name Null-terminated operation name to look up; if null the function
 * returns false.
 * @return true if the operation was found and printed; false if the name was
 * null or no matching operation exists.
 */
auto
print_operation_help(const char *name) -> bool
{
  if (name == nullptr) {
    return false;
  }

  auto match_result = find_operation(name);

  if (match_result.found) {
    std::cout << "\nOperation: "
              << to_underscore_convention(match_result.operation->name) << "\n"
              << "Category: " << match_result.operation->category << "\n"
              << "Description: " << match_result.operation->description << "\n";
    if (match_result.operation->requires_b) {
      std::cout << "Requires two geometries\n";
    }
    if (!match_result.operation->param_help.empty()) {
      std::cout << "\n" << match_result.operation->param_help << "\n";
    }
    return true;
  }

  std::cerr << "Unknown operation: " << name << "\n";
  return false;
}

/**
 * @brief Returns metadata for all registered operations.
 *
 * Produces a list of tuples describing every operation in the registry. Each
 * tuple contains:
 *  - name (std::string): operation identifier
 *  - category (std::string): grouping/category name
 *  - description (std::string): human-readable description
 *  - input (std::string): input specification (A, A,B, A,params)
 *  - output (std::string): output type (G, D, B, T)
 *
 * @return std::vector<std::tuple<std::string, std::string, std::string,
 * std::string, std::string>> Vector of operation metadata tuples in the order
 * they appear in the registry.
 */
auto
get_all_operations_info() -> std::vector<
    std::tuple<std::string, std::string, std::string, std::string, std::string>>
{
  std::vector<Operation> operations = OperationRegistry::all();

  std::vector<std::tuple<std::string, std::string, std::string, std::string,
                         std::string>>
      result;
  result.reserve(operations.size());

  for (const auto &operation : operations) {
    result.emplace_back(to_underscore_convention(operation.name),
                        operation.category, operation.description,
                        operation.input, operation.output);
  }

  return result;
}

/**
 * @brief Check if an operation requires a second geometry parameter.
 *
 * Looks up the operation by name in the operations registry and returns
 * whether it requires a second geometry parameter based on the requires_b
 * field in the operation definition.
 *
 * @param operation_name Name of the operation to check.
 * @return true if the operation requires a second geometry; false if it
 * doesn't require one or if the operation name is not found.
 */
auto
operation_requires_second_geometry(const std::string &operation_name) -> bool
{
  auto match_result = find_operation(operation_name);

  if (match_result.found) {
    return match_result.operation->requires_b;
  }

  return false; // Operation not found, assume it doesn't require second
                // geometry
}

} // namespace Operations

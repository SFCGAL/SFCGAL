// Copyright (c) 2026-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "operations_processing.hpp"
#include "SFCGAL/algorithm/Chamfer.h"
#include "SFCGAL/algorithm/eliminateSlivers.h"

#include "SFCGAL/Geometry.h"

namespace Operations {

auto
get_operations_processing() -> const std::vector<Operation> &
{
  static const std::vector<Operation> operations_processing = {
      {"chamfer", "Processing",
       "Apply a chamfer or fillet operation to a solid along an edge", true,
       "Parameters:\n"
       "  radius=VALUE: Primary radius/width (default: 0.1)\n"
       "  radius_y=VALUE: Secondary radius for asymmetric chamfer (default: "
       "-1, "
       "symmetric)\n"
       "  type=0|1: Operation type (0=flat, 1=round) (default: 0)\n"
       "  segments=INT: Number of segments for round fillet (default: 8)\n"
       "  epsilon=VALUE: Tolerance for edge matching (default: 1e-6)\n\n"
       "Examples:\n"
       "  sfcgalop -a \"SOLID((...))\" -b \"LINESTRING(...)\" chamfer "
       "\"radius=0.5\"\n"
       "  echo \"SOLID(...)\" | sfcgalop -b \"LINESTRING(...)\" chamfer "
       "\"type=1,radius=0.2,segments=16\"",
       "B", "",
       [](const std::string &args, const SFCGAL::Geometry *geom_a,
          const SFCGAL::Geometry *geom_b) -> std::optional<OperationResult> {
         if (!geom_a || !geom_b) {
           std::cerr
               << "chamfer error: requires two geometries (solid and edge)\n";
           return std::nullopt;
         }

         if (geom_b->geometryTypeId() != SFCGAL::TYPE_LINESTRING &&
             geom_b->geometryTypeId() != SFCGAL::TYPE_MULTILINESTRING) {
           std::cerr
               << "chamfer error: second geometry must be a LineString or "
                  "MultiLineString (edge)\n";
           return std::nullopt;
         }

         auto                              params = parse_params(args);
         SFCGAL::algorithm::ChamferOptions options;

         if (params.count("radius")) {
           options.radius = params["radius"];
         }
         if (params.count("radius_y")) {
           options.radius_y = params["radius_y"];
         }
         if (params.count("segments")) {
           options.segments = static_cast<int>(params["segments"]);
         }
         if (params.count("epsilon")) {
           options.epsilon = params["epsilon"];
         }
         if (params.count("type")) {
           int type_int = static_cast<int>(params["type"]);
           options.type = (type_int == 1)
                              ? SFCGAL::algorithm::ChamferType::ROUND
                              : SFCGAL::algorithm::ChamferType::FLAT;
         }

         try {
           return SFCGAL::algorithm::chamfer(*geom_a, *geom_b, options);
         } catch (const std::exception &e) {
           std::cerr << "chamfer error: " << e.what() << "\n";
           return std::nullopt;
         }
       }},

      {"eliminate_slivers", "Processing",
       "Collapse the sliver faces of a volumic geometry", false,
       "A boolean operation on coordinates that were rounded beforehand leaves "
       "needle faces as large as that rounding. They are valid, and degenerate "
       "as soon as the geometry is written with a coarser precision. This "
       "operation collapses the former and drops the latter.\n\n"
       "Parameters:\n"
       "  tolerance=VALUE: Longest edge a collapse may consume, in the unit of "
       "the coordinates\n"
       "                   (default: the rounding error of a double precision "
       "construction)\n\n"
       "Examples:\n"
       "  sfcgalop -a \"SOLID Z (...)\" eliminate_slivers\n"
       "  # a difference of geometries written with 6 decimals\n"
       "  sfcgalop -a \"SOLID Z (...)\" eliminate_slivers \"tolerance=1e-6\"",
       "A, params", "G",
       [](const std::string &args, const SFCGAL::Geometry *geom_a,
          const SFCGAL::Geometry *) -> std::optional<OperationResult> {
         auto         params = parse_params(args);
         const double tolerance =
             params.count("tolerance") != 0 ? params["tolerance"] : 0.0;
         return SFCGAL::algorithm::eliminateSlivers(*geom_a, tolerance);
       }}};

  return operations_processing;
}

} // namespace Operations

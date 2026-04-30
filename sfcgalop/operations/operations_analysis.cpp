// Copyright (c) 2026-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "operations_analysis.hpp"

#include "SFCGAL/Geometry.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/algorithm/orientation.h"
#include "SFCGAL/algorithm/partition_2.h"
#include "SFCGAL/algorithm/visibility.h"

#include <optional>

#include <CGAL/number_utils.h>

namespace Operations {

const std::vector<Operation> operations_analysis = {
    {"normal", "Analysis", "Compute surface normal vector for polygon/triangle",
     false, "", "A", "G",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       if (geom_a->geometryTypeId() == SFCGAL::TYPE_POLYGON ||
           geom_a->geometryTypeId() == SFCGAL::TYPE_TRIANGLE) {
         // Compute normal using first 3 points of polygon
         const auto &polygon =
             (geom_a->geometryTypeId() == SFCGAL::TYPE_POLYGON)
                 ? geom_a->as<SFCGAL::Polygon>()
                 : SFCGAL::Polygon(geom_a->as<SFCGAL::Triangle>());

         if (polygon.exteriorRing().numPoints() >= 3) {
           const auto &point0 = polygon.exteriorRing().pointN(0);
           const auto &point1 = polygon.exteriorRing().pointN(1);
           const auto &point2 = polygon.exteriorRing().pointN(2);

           // Calculate normal as cross product
           double v1x = CGAL::to_double(point1.x() - point0.x());
           double v1y = CGAL::to_double(point1.y() - point0.y());
           double v1z =
               point1.is3D() ? CGAL::to_double(point1.z() - point0.z()) : 0;

           double v2x = CGAL::to_double(point2.x() - point0.x());
           double v2y = CGAL::to_double(point2.y() - point0.y());
           double v2z =
               point2.is3D() ? CGAL::to_double(point2.z() - point0.z()) : 0;

           double nx = (v1y * v2z) - (v1z * v2y);
           double ny = (v1z * v2x) - (v1x * v2z);
           double nz = (v1x * v2y) - (v1y * v2x);

           return std::make_unique<SFCGAL::Point>(nx, ny, nz);
         }
       }
       return std::nullopt;
     }},

    {"orientation", "Analysis",
     "Determine polygon ring orientation (clockwise/counter-clockwise)", false,
     "", "A", "D",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       if (geom_a->geometryTypeId() == SFCGAL::TYPE_POLYGON) {
         const auto &polygon = geom_a->as<SFCGAL::Polygon>();
         // Check if polygon exterior ring is counter-clockwise
         bool ccw = SFCGAL::algorithm::isCounterClockWiseOriented(
             polygon.exteriorRing());
         return static_cast<double>(ccw ? 1 : -1);
       }
       return std::nullopt;
     }},

    {"partition", "Analysis", "Partition polygon into simpler pieces", false,
     "Parameters:\n  method=0|1|2|3 (default: 0)\n\nMethods:\n  0 = "
     "y_monotone: Creates y-monotone polygons\n  1 = approx_convex: "
     "Approximate convex partition\n  2 = greene_approx_convex: Greene's "
     "approximation algorithm\n  3 = optimal_convex: Optimal convex "
     "partition\n\nExample:\n  sfcgalop -a \"POLYGON((0 0,5 5,10 0,5 -2,0 "
     "0))\" partition \"method=1\"",
     "A, params", "G",
     [](const std::string &args, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto params = parse_params(args);
       int  method =
           static_cast<int>(params.count("method") ? params["method"] : 0);

       SFCGAL::algorithm::PartitionAlgorithm alg;
       switch (method) {
       case 1:
         alg = SFCGAL::algorithm::approx_convex;
         break;
       case 2:
         alg = SFCGAL::algorithm::greene_approx_convex;
         break;
       case 3:
         alg = SFCGAL::algorithm::optimal_convex;
         break;
       default:
         alg = SFCGAL::algorithm::y_monotone;
       }

       return SFCGAL::algorithm::partition_2(*geom_a, alg);
     }},

    {"visibility", "Analysis",
     "Compute visibility polygon from a point in polygon", false,
     "Parameters:\n  x=X_COORD: X coordinate of viewpoint\n  y=Y_COORD: Y "
     "coordinate of viewpoint\n\nExample:\n  sfcgalop -a \"POLYGON((0 0,10 "
     "0,10 10,0 10,0 0))\" visibility \"x=5,y=5\"",
     "A, params", "G",
     [](const std::string &args, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       if (geom_a->geometryTypeId() == SFCGAL::TYPE_POLYGON) {
         const auto   &polygon = geom_a->as<SFCGAL::Polygon>();
         auto          params  = parse_params(args);
         SFCGAL::Point point(params["x"], params["y"]);
         return SFCGAL::algorithm::visibility(polygon, point);
       }
       return std::nullopt;
     }}};

} // namespace Operations

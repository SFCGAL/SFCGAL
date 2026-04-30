// Copyright (c) 2026-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "operations_metrics.hpp"

#include "SFCGAL/algorithm/area.h"
#include "SFCGAL/algorithm/distance.h"
#include "SFCGAL/algorithm/distance3d.h"
#include "SFCGAL/algorithm/length.h"
#include "SFCGAL/algorithm/volume.h"

namespace Operations {

const std::vector<Operation> operations_metrics = {
    {"area", "Metrics", "Calculate the 2D area of a geometry", false,
     "No parameters required.\n\nExample:\n  sfcgalop -a \"POLYGON((0 0,3 0,3 "
     "4,0 4,0 0))\" area",
     "A", "D",
     [](const std::string &, const SFCGAL::Geometry *geom,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       return SFCGAL::algorithm::area(*geom);
     }},

    {"area_3d", "Metrics", "Calculate the 3D surface area of a geometry", false,
     "No parameters required.\n\nExample:\n  sfcgalop -a \"POLYGON Z((0 0 0,3 "
     "0 0,3 4 2,0 4 2,0 0 0))\" area_3d",
     "A", "D",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       return SFCGAL::algorithm::area3D(*geom_a);
     }},

    {"distance", "Metrics",
     "Calculate the 2D minimum distance between two geometries", true,
     "No parameters required.\nRequires two geometries.\n\nExample:\n  "
     "sfcgalop -a \"POINT(0 0)\" -b \"POINT(3 4)\" distance",
     "A, B", "D",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *geom_b) -> std::optional<OperationResult> {
       if (!geom_b) {
         return std::nullopt;
       }
       return SFCGAL::algorithm::distance(*geom_a, *geom_b);
     }},

    {"distance_3d", "Metrics",
     "Calculate the 3D minimum distance between two geometries", true,
     "No parameters required.\nRequires two geometries.\n\nExample:\n  "
     "sfcgalop -a \"POINT Z(0 0 0)\" -b \"POINT Z(3 4 5)\" distance_3d",
     "A, B", "D",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *geom_b) -> std::optional<OperationResult> {
       if (!geom_b) {
         return std::nullopt;
       }
       return SFCGAL::algorithm::distance3D(*geom_a, *geom_b);
     }},

    {"length", "Metrics", "Calculate the 2D length of linear geometries", false,
     "No parameters required.\n\nExample:\n  sfcgalop -a \"LINESTRING(0 0,3 "
     "4,6 0)\" length",
     "A", "D",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       return SFCGAL::algorithm::length(*geom_a);
     }},

    {"length_3d", "Metrics", "Calculate the 3D length of linear geometries",
     false,
     "No parameters required.\n\nExample:\n  sfcgalop -a \"LINESTRING Z(0 0 "
     "0,3 4 2,6 0 1)\" length_3d",
     "A", "D",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       return SFCGAL::algorithm::length3D(*geom_a);
     }},

    {"volume", "Metrics", "Calculate the 3D volume of a solid geometry", false,
     "No parameters required.\nOnly works with solid geometries (SOLID, "
     "POLYHEDRALSURFACE).\n\nExample:\n  sfcgalop -a \"SOLID(...)\" volume",
     "A", "D",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       return CGAL::to_double(SFCGAL::algorithm::volume(*geom_a));
     }}};

} // namespace Operations

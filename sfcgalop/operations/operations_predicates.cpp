// Copyright (c) 2026-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "operations_predicates.hpp"

#include "SFCGAL/algorithm/covers.h"
#include "SFCGAL/algorithm/intersects.h"
#include "SFCGAL/algorithm/isClosed.h"
#include "SFCGAL/algorithm/isSimple.h"
#include "SFCGAL/algorithm/isValid.h"

namespace Operations {

const std::vector<Operation> operations_predicates = {
    {"covers", "Predicates", "Test if geometry A completely covers geometry B",
     true, "", "A, B", "B",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *geom_b) -> std::optional<OperationResult> {
       if (!geom_b) {
         return std::nullopt;
       }
       return SFCGAL::algorithm::covers(*geom_a, *geom_b);
     }},

    {"intersects", "Predicates", "Test if two geometries intersect in 2D", true,
     "", "A, B", "B",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *geom_b) -> std::optional<OperationResult> {
       if (!geom_b) {
         return std::nullopt;
       }
       return SFCGAL::algorithm::intersects(*geom_a, *geom_b);
     }},

    {"intersects_3d", "Predicates", "Test if two geometries intersect in 3D",
     true, "", "A, B", "B",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *geom_b) -> std::optional<OperationResult> {
       if (!geom_b) {
         return std::nullopt;
       }
       return SFCGAL::algorithm::intersects3D(*geom_a, *geom_b);
     }},

    {"is_3d", "Predicates", "Test if geometry has Z coordinates", false, "",
     "A", "B",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       return geom_a->is3D();
     }},

    {"is_closed", "Predicates", "Test if linear geometry forms a closed ring",
     false, "", "A", "B",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       return static_cast<bool>(SFCGAL::algorithm::isClosed(*geom_a));
     }},

    {"is_empty", "Predicates", "Test if geometry contains no points", false, "",
     "A", "B",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       return geom_a->isEmpty();
     }},

    {"is_measured", "Predicates",
     "Test if geometry has measure (M) coordinates", false, "", "A", "B",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       return geom_a->isMeasured();
     }},

    {"is_simple", "Predicates", "Test if geometry has no self-intersections",
     false, "", "A", "B",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       return static_cast<bool>(SFCGAL::algorithm::isSimple(*geom_a));
     }},

    {"is_valid", "Predicates", "Test if geometry is topologically valid", false,
     "", "A", "B",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       return static_cast<bool>(SFCGAL::algorithm::isValid(*geom_a));
     }}};

} // namespace Operations

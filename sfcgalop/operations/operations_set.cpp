// Copyright (c) 2026-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "operations_set.hpp"

#include "SFCGAL/algorithm/difference.h"
#include "SFCGAL/algorithm/intersection.h"
#include "SFCGAL/algorithm/length.h"
#include "SFCGAL/algorithm/union.h"

namespace Operations {

const std::vector<Operation> operations_set = {
    {"difference", "Set Operations", "Compute geometry A minus geometry B",
     true, "", "A, B", "G",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *geom_b) -> std::optional<OperationResult> {
       if (!geom_b) {
         return std::nullopt;
       }
       return SFCGAL::algorithm::difference(*geom_a, *geom_b);
     }},

    {"difference_3d", "Set Operations",
     "Compute 3D geometry A minus geometry B", true, "", "A, B", "G",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *geom_b) -> std::optional<OperationResult> {
       if (!geom_b) {
         return std::nullopt;
       }
       return SFCGAL::algorithm::difference3D(*geom_a, *geom_b);
     }},

    {"intersection", "Set Operations",
     "Compute the geometric intersection of two geometries", true, "", "A, B",
     "G",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *geom_b) -> std::optional<OperationResult> {
       if (!geom_b) {
         return std::nullopt;
       }
       return SFCGAL::algorithm::intersection(*geom_a, *geom_b);
     }},

    {"intersection_3d", "Set Operations",
     "Compute the 3D geometric intersection of two geometries", true, "",
     "A, B", "G",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *geom_b) -> std::optional<OperationResult> {
       if (!geom_b) {
         return std::nullopt;
       }
       return SFCGAL::algorithm::intersection3D(*geom_a, *geom_b);
     }},

    {"union", "Set Operations", "Compute the geometric union of two geometries",
     true, "", "A, B", "G",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *geom_b) -> std::optional<OperationResult> {
       if (!geom_b) {
         return std::nullopt;
       }
       return SFCGAL::algorithm::union_(*geom_a, *geom_b);
     }},

    {"union_3d", "Set Operations",
     "Compute the 3D geometric union of two geometries", true, "", "A, B", "G",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *geom_b) -> std::optional<OperationResult> {
       if (!geom_b) {
         return std::nullopt;
       }
       return SFCGAL::algorithm::union3D(*geom_a, *geom_b);
     }}};

} // namespace Operations

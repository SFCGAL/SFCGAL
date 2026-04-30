// Copyright (c) 2026-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "operations_conversions.hpp"
#include "../constructors.hpp"

namespace Operations {

const std::vector<Operation> operations_conversions = {
    {"to_solid", "Conversions", "Convert a PolyhedralSurface to a Solid", false,
     "", "A", "",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       // Clone the input geometry to pass ownership to make_solid
       auto geom_copy = geom_a->clone();
       return Constructors::make_solid(std::move(geom_copy));
     }}};

} // namespace Operations

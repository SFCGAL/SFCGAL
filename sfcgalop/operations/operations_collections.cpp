// Copyright (c) 2026-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "operations_collections.hpp"

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/algorithm/collectionExtract.h"
#include "SFCGAL/algorithm/collectionHomogenize.h"
#include "SFCGAL/algorithm/collectionToMulti.h"

#include <memory>

namespace Operations {

const std::vector<Operation> operations_collections = {
    {"collect", "Collections", "Combine two geometries into a collection", true,
     "", "A, B", "G",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *geom_b) -> std::optional<OperationResult> {
       if (!geom_b) {
         return std::nullopt;
       }
       auto collection = std::make_unique<SFCGAL::GeometryCollection>();
       collection->addGeometry(*geom_a);
       collection->addGeometry(*geom_b);
       return collection;
     }},

    {"collection_extract", "Collections",
     "Extract polygons from geometry collection", false, "", "A", "G",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto copy = geom_a->clone();
       return SFCGAL::algorithm::collectionExtractPolygons(std::move(copy));
     }},

    {"collection_homogenize", "Collections",
     "Convert collection to appropriate multi-type", false, "", "A", "G",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto copy = geom_a->clone();
       return SFCGAL::algorithm::collectionHomogenize(std::move(copy));
     }},

    {"collection_to_multi", "Collections",
     "Convert collection to multi-geometry type", false, "", "A", "G",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto copy = geom_a->clone();
       return SFCGAL::algorithm::collectionToMulti(std::move(copy));
     }}};

} // namespace Operations

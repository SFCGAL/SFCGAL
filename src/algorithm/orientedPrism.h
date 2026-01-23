// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_ORIENTED_PRISM_H_
#define SFCGAL_ALGORITHM_ORIENTED_PRISM_H_

#include "SFCGAL/Geometry.h"
#include "SFCGAL/config.h"

#include <memory>

namespace SFCGAL::algorithm {

/**
 * @brief Constructs a prism (or series of prisms) oriented along a linestring
 * using incident faces of a solid to determine the base orientation.
 *
 * @param solid The solid containing the edges
 * @param path The linestring path (edges of the solid)
 * @param distance1 Offset distance on the first incident face
 * @param distance2 Offset distance on the second incident face
 * @return A unique_ptr to the generated geometry (Solid or MultiSolid)
 */
SFCGAL_API auto
orientedPrism(const Geometry &solid, const Geometry &path, double distance1,
              double distance2) -> std::unique_ptr<Geometry>;

} // namespace SFCGAL::algorithm

#endif // SFCGAL_ALGORITHM_ORIENTED_PRISM_H_

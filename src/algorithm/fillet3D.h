// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_FILLET3D_H_
#define SFCGAL_ALGORITHM_FILLET3D_H_

#include "SFCGAL/Geometry.h"
#include "SFCGAL/config.h"

#include <memory>

namespace SFCGAL::algorithm {

/**
 * @brief Computes the 3D fillet of a solid along specified edges (by radius).
 *
 * A fillet operation creates rounded transitions along edges of a solid by
 * generating a "rounded scraper" volume (prism minus sphere) and subtracting
 * it from the solid. The sphere radius is r, and symmetric distances are
 * computed as distance = r / sqrt(2).
 *
 * @param solid The base solid geometry to fillet (Solid or PolyhedralSurface)
 * @param edges The edges along which to apply the fillet (LineString or
 * MultiLineString)
 * @param radius The radius of the fillet (sphere radius = sqrt(d1² + d2²))
 * @param num_subdivisions Number of sphere subdivisions for rounding quality
 * (default=2)
 * @return The filleted solid as a unique_ptr<Geometry>
 *
 * @throws SFCGAL::Exception if the input geometry types are invalid
 */
SFCGAL_API auto
fillet3D(const Geometry &solid, const Geometry &edges, double radius,
         unsigned int num_subdivisions = 2) -> std::unique_ptr<Geometry>;

/**
 * @brief Computes the 3D fillet of a solid along specified edges (by
 * distances).
 *
 * A fillet operation creates rounded transitions along edges of a solid by
 * generating a "rounded scraper" volume (prism minus sphere) and subtracting
 * it from the solid.
 *
 * @param solid The base solid geometry to fillet (Solid or PolyhedralSurface)
 * @param edges The edges along which to apply the fillet (LineString or
 * MultiLineString)
 * @param distance1 The distance to offset along the first incident face
 * @param distance2 The distance to offset along the second incident face
 * @param num_subdivisions Number of sphere subdivisions for rounding quality
 * (default=2)
 * @return The filleted solid as a unique_ptr<Geometry>
 *
 * @throws SFCGAL::Exception if the input geometry types are invalid
 */
SFCGAL_API auto
fillet3D_asymmetric(const Geometry &solid, const Geometry &edges,
                    double distance1, double distance2,
                    unsigned int num_subdivisions = 2)
    -> std::unique_ptr<Geometry>;

} // namespace SFCGAL::algorithm

#endif // SFCGAL_ALGORITHM_FILLET3D_H_

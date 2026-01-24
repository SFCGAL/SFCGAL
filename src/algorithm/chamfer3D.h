// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_CHAMFER3D_H_
#define SFCGAL_ALGORITHM_CHAMFER3D_H_

#include "SFCGAL/Geometry.h"
#include "SFCGAL/config.h"

#include <memory>

namespace SFCGAL::algorithm {

/**
 * @brief Computes the 3D chamfer of a solid along specified edges.
 *
 * A chamfer operation removes triangular prisms along edges of a solid to
 * create beveled transitions between incident faces.
 *
 * @param solid The base solid geometry to chamfer (Solid or PolyhedralSurface)
 * @param edges The edges along which to apply the chamfer (LineString or
 * MultiLineString)
 * @param distance1 The distance to offset along the first incident face
 * @param distance2 The distance to offset along the second incident face
 * @return The chamfered solid as a unique_ptr<Geometry>
 *
 * @throws SFCGAL::Exception if the input geometry types are invalid
 */
SFCGAL_API auto
chamfer3D(const Geometry &solid, const Geometry &edges, double distance1,
          double distance2) -> std::unique_ptr<Geometry>;

/**
 * @brief Computes the symmetric 3D chamfer of a solid along specified edges.
 *
 * @param solid The base solid geometry
 * @param edges The edges along which to apply the chamfer
 * @param distance The symmetric offset distance for both incident faces
 * @return The chamfered solid
 */
inline auto
chamfer3D(const Geometry &solid, const Geometry &edges, double distance)
    -> std::unique_ptr<Geometry>
{
  return chamfer3D(solid, edges, distance, distance);
}

} // namespace SFCGAL::algorithm

#endif // SFCGAL_ALGORITHM_CHAMFER3D_H_

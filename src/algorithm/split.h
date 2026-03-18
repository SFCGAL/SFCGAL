// Copyright (c) 2026-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_SPLIT_H_
#define SFCGAL_ALGORITHM_SPLIT_H_

#include "SFCGAL/Geometry.h"
#include "SFCGAL/Kernel.h"
#include "SFCGAL/Point.h"

namespace SFCGAL::algorithm {

/**
 * Split the given geometry with a plane defined by a point and a normal vector.
 * Based on CGAL::PMP::split operator
 *
 * Behavior by geometry type:
 * - POINT, LINESTRING, NURBSCURVE, MULTIPOINT, MULTILINESTRING, TRIANGLE: not
 * handled, throws an exception
 * - POLYGON, TRIANGULATEDSURFACE, POLYHEDRALSURFACE: splits the geometry
 * - MULTI types and GEOMETRYCOLLECTION: iterates over the geometries and splits
     each geometry
 *
 * @param geometry The input geometry to split.
 * @param planePoint A point on the splitting plane.
 * @param planeNormal The normal vector defining the plane orientation.
 * @param closeGeometries If true, ensures resulting geometries are closed.
 *   In the case of a SOLID, this parameter is always forced to true,
 *   since a solid is a closed geometry.
 *
 * @return A unique pointer to a GeometryCollection containing the split parts
 * or a unique pointer to an empty GeometryCollection if the plane does not
 * intersect the geometry.
 *
 * @since 2.3
 */
SFCGAL_API auto
split(const Geometry &geometry, const Point &planePoint,
      const Kernel::Vector_3 &planeNormal, bool closeGeometries = false)
    -> std::unique_ptr<GeometryCollection>;

} // namespace SFCGAL::algorithm

#endif // SFCGAL_ALGORITHM_SPLIT_H_

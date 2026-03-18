// Copyright (c) 2026-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_SPLIT_3D_H_
#define SFCGAL_ALGORITHM_SPLIT_3D_H_

#include "SFCGAL/Geometry.h"
#include "SFCGAL/Kernel.h"
#include "SFCGAL/Point.h"

namespace SFCGAL::algorithm {

/**
 * Split the given geometry with a plane defined by a point and a normal vector.
 * Based on CGAL::PMP::split operator
 * https://doc.cgal.org/latest/Polygon_mesh_processing/group__PMP__corefinement__grp.html
 *
 * Only the TIN, PolyhedralSurface, and Solid types are currently supported.
 *
 * @param geometry The input geometry to split.
 * @param planePoint A point on the splitting plane.
 * @param planeNormal The normal vector defining the plane orientation.
 * @param closeGeometries If true, ensures resulting geometries are closed.
 *   In the case of a SOLID, this parameter is always forced to true,
 *   since a solid is a closed geometry.
 *
 * @return A unique pointer to a GeometryCollection of the split parts
 * or the original geometry if the plane does not intersect it.
 *
 * @since 2.3
 */
SFCGAL_API auto
split3D(const Geometry &geometry, const Point &planePoint,
        const Kernel::Vector_3 &planeNormal, bool closeGeometries = false)
    -> std::unique_ptr<GeometryCollection>;

} // namespace SFCGAL::algorithm

#endif // SFCGAL_ALGORITHM_SPLIT_3D_H_

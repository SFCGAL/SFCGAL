// Copyright (c) 2026-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_DETAIL_ALGORITHM_SPLIT_3D_H_
#define SFCGAL_DETAIL_ALGORITHM_SPLIT_3D_H_

#include <concepts>
#include <memory>

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/Kernel.h"

namespace SFCGAL::algorithm::detail {

template <typename T>
concept SurfaceMeshable = requires(const T &geometry) {
  { geometry.toSurfaceMesh() } -> std::same_as<Surface_mesh_3>;
};

/**
 * Split the given geometry with a plane.
 * Based on CGAL::PMP::split operator
 *
 * @param geometry The input geometry to split.
 * @param plane The splitting plane
 * @param closeGeometries If true, ensures resulting geometries are closed.
 *   In the case of a SOLID, this parameter is always forced to true,
 *   since a solid is a closed geometry.
 *
 * @return A unique pointer to a GeometryCollection of the split parts
 * or the original geometry if the plane does not intersect it.
 */
template <SurfaceMeshable geomType>
auto
split3D(const geomType &geometry, const CGAL::Plane_3<Kernel> &plane,
        bool closeGeometries) -> std::unique_ptr<GeometryCollection>;

} // namespace SFCGAL::algorithm::detail

#endif // SFCGAL_DETAIL_ALGORITHM_SPLIT_3D_H_

// Copyright (c) 2026-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_DETAIL_SPLIT_H_
#define SFCGAL_ALGORITHM_DETAIL_SPLIT_H_

#include <memory>

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/Kernel.h"

namespace SFCGAL::algorithm::detail {

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
 * @return A unique pointer to a GeometryCollection containing the split parts
 * or a unique pointer to an empty GeometryCollection if the plane does not
 * intersect the geometry.
 */
template <typename geomType>
auto
split(const geomType &geometry, const CGAL::Plane_3<Kernel> &plane,
      bool closeGeometries) -> std::unique_ptr<GeometryCollection>;

} // namespace SFCGAL::algorithm::detail

#endif // SFCGAL_ALGORITHM_DETAIL_SPLIT_H_

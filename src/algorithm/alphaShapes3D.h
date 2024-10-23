// Copyright (c) 2026-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_ALPHASHAPES3D_H_
#define SFCGAL_ALGORITHM_ALPHASHAPES3D_H_

#include "SFCGAL/Geometry.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/config.h"

namespace SFCGAL::algorithm {

/**
 * Compute the 3D alpha shape for a geometry with optimal alpha value
 * https://doc.cgal.org/latest/Alpha_shapes_3/index.html
 * @ingroup public_api
 * @since 2.3
 * @param geometry input geometry
 * @param allowCavities Whether to allow cavities in the result (default: false)
 * @return A PolyhedralSurface representing the optimal alpha shape
 */
SFCGAL_API auto
alphaShapes3D(const Geometry &geometry, bool allowCavities = false)
    -> std::unique_ptr<PolyhedralSurface>;

} // namespace SFCGAL::algorithm

#endif // SFCGAL_ALGORITHM_ALPHASHAPES3D_H_

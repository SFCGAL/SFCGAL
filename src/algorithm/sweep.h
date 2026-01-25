// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_SWEEP_H_
#define SFCGAL_ALGORITHM_SWEEP_H_

#include "SFCGAL/Geometry.h"
#include "SFCGAL/config.h"

#include <memory>

namespace SFCGAL::algorithm {

/**
 * @brief Sweeps a 2D cross-section geometry along a path (LineString)
 *
 * This creates a 3D surface by moving a 2D profile along a path, similar to
 * extrusion but along an arbitrary path with proper corner handling.
 *
 * The cross-section is positioned perpendicular to the path at each point,
 * and bisector planes are used to join consecutive segments smoothly.
 *
 * @param path The path along which to sweep (LineString or MultiLineString)
 * @param cross_section The 2D profile to sweep (Point, LineString, or Polygon)
 * @param close_ends Whether to close the ends with caps (default: true)
 * @return A PolyhedralSurface representing the swept volume
 *
 * @note The cross_section should be defined in a local 2D coordinate system
 *       and will be oriented perpendicular to the path automatically
 *
 * Examples:
 * - Sweep a circle along a path → pipe/tube
 * - Sweep a triangle along an edge → chamfer volume
 * - Sweep a "triangle minus arc" → fillet volume
 * - Sweep a rectangle → ribbon/band
 *
 * @throws SFCGAL::Exception if inputs are invalid
 */
SFCGAL_API auto
sweep(const Geometry &path, const Geometry &cross_section,
      bool close_ends = true) -> std::unique_ptr<Geometry>;

} // namespace SFCGAL::algorithm

#endif // SFCGAL_ALGORITHM_SWEEP_H_

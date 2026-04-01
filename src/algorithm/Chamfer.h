// Copyright (c) 2025-2026, SFCGAL Team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_CHAMFER_H_
#define SFCGAL_ALGORITHM_CHAMFER_H_

#include <SFCGAL/Geometry.h>
#include <SFCGAL/config.h>

#include <memory>
#include <vector>

namespace SFCGAL::algorithm {

/**
 * @brief Type of chamfer operation
 *
 * A chamfer removes material along a sharp edge of a solid.
 * - **FLAT**: Replaces the edge with a single planar face (a bevel).
 * - **ROUND**: Replaces the edge with a smooth curved surface (a fillet),
 *   approximated by a series of planar facets.
 */
enum class ChamferType {
  FLAT, ///< Flat chamfer: single planar bevel
  ROUND ///< Rounded fillet: circular-arc surface (approximated)
};

/**
 * @brief Options for chamfer/fillet operations
 */
struct ChamferOptions {
  /// Type of operation: FLAT (chamfer) or ROUND (fillet)
  ChamferType type = ChamferType::FLAT;

  /// Distance from the edge along each face surface.
  /// For FLAT: length of the chamfer leg on the first face.
  /// For ROUND: radius of the circular arc tangent to both faces.
  /// Must be positive.
  double radius = 0.1;

  /// Secondary distance for asymmetric FLAT chamfers: leg length on the
  /// second face. If negative, defaults to @c radius (symmetric chamfer).
  /// Ignored for ROUND type.
  double radius_y = -1.0;

  /// Number of linear segments approximating the fillet arc (ROUND only).
  /// Typical values: 4 (coarse), 8 (default), 16+ (smooth).
  /// Must be >= 1. Ignored for FLAT type.
  int segments = 8;

  /// Tolerance (in coordinate units) for matching input edge endpoints to
  /// the solid's mesh vertices. Default 1e-3 accommodates numeric drift
  /// from coordinate transformations.
  double epsilon = 1e-3;
};

/**
 * @brief Apply a chamfer or fillet to a solid along specified edges.
 *
 * A **chamfer** (FLAT) replaces a sharp edge with a flat bevel. A **fillet**
 * (ROUND) replaces it with a smooth circular-arc surface. The cutter profile
 * automatically adapts to the actual dihedral angle between the two faces
 * meeting at the edge, supporting arbitrary convex edges (not just 90°).
 *
 * **How it works:** For each edge, the function identifies the two incident
 * faces, measures their dihedral angle, builds a cutting profile, sweeps it
 * along the edge, and subtracts the result from the solid.
 *
 * **Multi-edge:** When a MultiLineString or multi-segment LineString is given,
 * each segment is processed independently (each mesh edge has its own dihedral
 * angle). The individual cutters are unioned before a single boolean difference.
 *
 * **Silent skipping:** Edges that fail validation (not found in mesh, concave,
 * angle out of range) are skipped with a warning log, not thrown as exceptions.
 * If all edges are skipped, the original solid is returned unmodified.
 *
 * @par Limitations
 * - Opening angle must be between ~10° and ~170°. Edges outside this range
 *   are skipped.
 * - Reflex (concave) edges are skipped.
 * - Both incident faces must be planar.
 * - The radius must be small enough that the cutter fits within the faces.
 *
 * @param solid Input 3D Solid or PolyhedralSurface.
 * @param edge LineString or MultiLineString coinciding with solid edges.
 * @param options Chamfer configuration (radius, type, segments, etc.)
 * @return Modified solid, or clone of input if no edges were chamfered.
 *
 * @throws std::invalid_argument If @p solid is not a Solid/PolyhedralSurface.
 * @throws std::invalid_argument If @p edge is not a LineString/MultiLineString.
 * @throws std::invalid_argument If options.radius <= 0.
 * @throws std::invalid_argument If options.type is ROUND and segments < 1.
 *
 * @note Per-edge errors (not found, concave, angle out of range) are logged
 *   as warnings via SFCGAL_WARNING and do not throw.
 */
SFCGAL_API auto
chamfer(const Geometry &solid, const Geometry &edge,
        const ChamferOptions &options = {}) -> std::unique_ptr<Geometry>;

} // namespace SFCGAL::algorithm

#endif // SFCGAL_ALGORITHM_CHAMFER_H_

// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/difference.h"

#include "SFCGAL/Exception.h"
#include "SFCGAL/Geometry.h"
#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"

#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>

namespace SFCGAL::algorithm {

namespace {

using Nef_polyhedron_3 = CGAL::Nef_polyhedron_3<Kernel>;
using Polyhedron_3     = CGAL::Polyhedron_3<Kernel>;

/**
 * @brief Convert SFCGAL::Geometry to Nef_polyhedron_3
 *
 * Handles all geometry types by converting to Polyhedron_3 first.
 * This reuses existing conversion logic in SFCGAL.
 *
 * @param g The geometry to convert
 * @return Nef_polyhedron_3 representation
 */
auto
geometryToNef(const Geometry &g) -> Nef_polyhedron_3
{
  if (g.isEmpty()) {
    return {};
  }

  // For 3D volumes (Solid, PolyhedralSurface), convert via Polyhedron_3
  if (g.geometryTypeId() == TYPE_SOLID) {
    const auto &solid   = g.as<Solid>();
    auto        polyPtr = solid.exteriorShell().toPolyhedron_3<Polyhedron_3>();
    if (polyPtr && !polyPtr->is_empty()) {
      Nef_polyhedron_3 result(*polyPtr);

      // Subtract interior shells (voids)
      for (size_t i = 0; i < solid.numInteriorShells(); ++i) {
        auto interiorPtr =
            solid.interiorShellN(i).toPolyhedron_3<Polyhedron_3>();
        if (interiorPtr && !interiorPtr->is_empty()) {
          result -= Nef_polyhedron_3(*interiorPtr);
        }
      }
      return result;
    }
  } else if (g.geometryTypeId() == TYPE_POLYHEDRALSURFACE) {
    auto polyPtr = g.as<PolyhedralSurface>().toPolyhedron_3<Polyhedron_3>();
    if (polyPtr && !polyPtr->is_empty()) {
      return Nef_polyhedron_3(*polyPtr);
    }
  } else if (g.is<GeometryCollection>()) {
    // Union all components
    const auto      &coll = g.as<GeometryCollection>();
    Nef_polyhedron_3 result;
    for (size_t i = 0; i < coll.numGeometries(); ++i) {
      Nef_polyhedron_3 nef = geometryToNef(coll.geometryN(i));
      if (!nef.is_empty()) {
        if (result.is_empty()) {
          result = nef;
        } else {
          result += nef;
        }
      }
    }
    return result;
  }

  return {};
}

/**
 * @brief Convert Nef_polyhedron_3 back to SFCGAL::Geometry
 *
 * @param nef The Nef polyhedron to convert
 * @return Geometry as PolyhedralSurface or empty GeometryCollection
 */
auto
nefToGeometry(const Nef_polyhedron_3 &nef) -> std::unique_ptr<Geometry>
{
  if (nef.is_empty()) {
    return std::make_unique<GeometryCollection>();
  }

  Polyhedron_3 poly;
  nef.convert_to_polyhedron(poly);

  if (poly.is_empty()) {
    return std::make_unique<GeometryCollection>();
  }

  return std::make_unique<PolyhedralSurface>(poly);
}

} // anonymous namespace

/**
 * @brief 3D Boolean difference using Nef polyhedra
 *
 * This function performs 3D difference operations using CGAL's
 * Nef_polyhedron_3, which provides robust boolean operations with exact
 * arithmetic.
 *
 * Advantages over the MarkedPolyhedron/corefine approach:
 * - Exact arithmetic throughout the operation
 * - No vertex coherence issues
 * - Handles all degenerate cases robustly
 *
 * @param geometry1 The first geometry (A)
 * @param geometry2 The second geometry (B)
 * @return A - B as a unique_ptr<Geometry>
 */
auto
difference3D_nef(const Geometry &geometry1, const Geometry &geometry2)
    -> std::unique_ptr<Geometry>
{
  Nef_polyhedron_3 nefA = geometryToNef(geometry1);
  Nef_polyhedron_3 nefB = geometryToNef(geometry2);

  if (nefA.is_empty()) {
    return std::make_unique<GeometryCollection>();
  }

  if (nefB.is_empty()) {
    return geometry1.clone();
  }

  // Perform boolean difference: A - B
  Nef_polyhedron_3 result = nefA - nefB;

  return nefToGeometry(result);
}

} // namespace SFCGAL::algorithm

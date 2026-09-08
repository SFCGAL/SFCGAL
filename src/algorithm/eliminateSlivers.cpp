// Copyright (c) 2026-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/eliminateSlivers.h"

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiSolid.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/TriangulatedSurface.h"
#include "SFCGAL/detail/algorithm/collapseSlivers.h"

namespace SFCGAL::algorithm {

namespace {

/// Drops the points a ring repeats in a row.
auto
withoutRepeatedPoints(const LineString &ring) -> LineString
{
  LineString result;
  for (std::size_t i = 0; i < ring.numPoints(); ++i) {
    if (result.isEmpty() || !(result.endPoint() == ring.pointN(i))) {
      result.addPoint(ring.pointN(i));
    }
  }
  return result;
}

/**
 * Drops the patches a sliver already collapsed into, which a coarse write
 * leaves behind: writing both ends of a needle to the same coordinates folds
 * its ring onto itself. Such a patch has no plane, so it has to go before the
 * shell can be triangulated at all. It carried no area, and the two halves of
 * the edge it kept cancel each other, so the shell stays closed.
 */
auto
withoutFoldedPatches(const PolyhedralSurface &shell) -> PolyhedralSurface
{
  PolyhedralSurface result;

  for (std::size_t p = 0; p < shell.numPatches(); ++p) {
    const Polygon &patch = shell.patchN(p);

    LineString exteriorRing = withoutRepeatedPoints(patch.exteriorRing());
    if (exteriorRing.numPoints() < 4) {
      continue;
    }

    Polygon repaired(exteriorRing);
    for (std::size_t r = 0; r < patch.numInteriorRings(); ++r) {
      LineString interiorRing = withoutRepeatedPoints(patch.interiorRingN(r));
      if (interiorRing.numPoints() >= 4) {
        repaired.addInteriorRing(interiorRing);
      }
    }
    result.addPatch(repaired);
  }
  return result;
}

/**
 * Rebuilds a shell through a MarkedPolyhedron, which is what the collapse
 * works on. The round trip triangulates, which is what a boolean operation
 * hands out anyway.
 */
auto
repairShell(const PolyhedralSurface &shell, double maxEdgeLength)
    -> PolyhedralSurface
{
  const PolyhedralSurface unfolded = withoutFoldedPatches(shell);

  std::unique_ptr<detail::MarkedPolyhedron> mesh =
      unfolded.toPolyhedron_3<detail::MarkedPolyhedron>();
  detail::algorithm::collapseSlivers(*mesh, maxEdgeLength);
  return {*mesh, false};
}

} // namespace

auto
eliminateSlivers(const Geometry &geometry, double maxEdgeLength)
    -> std::unique_ptr<Geometry>
{
  switch (geometry.geometryTypeId()) {
  case TYPE_TRIANGULATEDSURFACE:
  case TYPE_POLYHEDRALSURFACE: {
    const PolyhedralSurface shell(geometry.clone());
    return std::make_unique<PolyhedralSurface>(
        repairShell(shell, maxEdgeLength));
  }

  case TYPE_SOLID: {
    const Solid &solid  = geometry.as<Solid>();
    auto         result = std::make_unique<Solid>(
        repairShell(solid.exteriorShell(), maxEdgeLength));
    for (std::size_t i = 0; i < solid.numInteriorShells(); ++i) {
      result->addInteriorShell(
          repairShell(solid.interiorShellN(i), maxEdgeLength));
    }
    return result;
  }

  case TYPE_MULTISOLID:
  case TYPE_GEOMETRYCOLLECTION: {
    const auto &collection = geometry.as<GeometryCollection>();
    auto        result     = std::make_unique<GeometryCollection>();
    for (std::size_t i = 0; i < collection.numGeometries(); ++i) {
      result->addGeometry(
          eliminateSlivers(collection.geometryN(i), maxEdgeLength));
    }
    return result;
  }

  default:
    return std::unique_ptr<Geometry>(geometry.clone());
  }
}

} // namespace SFCGAL::algorithm

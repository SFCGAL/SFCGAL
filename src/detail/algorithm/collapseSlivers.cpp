// Copyright (c) 2026-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/detail/algorithm/collapseSlivers.h"

#include "SFCGAL/detail/tools/Log.h"

#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/repair_degeneracies.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>

#include <cmath>
#include <limits>
#include <sstream>

namespace PMP = CGAL::Polygon_mesh_processing;

namespace SFCGAL::detail::algorithm {

namespace {

/// Ratio of the longest edge over the shortest one above which CGAL calls a
/// face a needle, and its own default.
constexpr double NEEDLE_RATIO = 4.0;

auto
boundingBoxDiagonal(const MarkedPolyhedron &mesh) -> double
{
  const CGAL::Bbox_3 box = PMP::bbox(mesh);
  return std::sqrt(CGAL::square(box.xmax() - box.xmin()) +
                   CGAL::square(box.ymax() - box.ymin()) +
                   CGAL::square(box.zmax() - box.zmin()));
}

/**
 * Reports whether the mesh holds a face the repair below is meant to remove:
 * a needle whose shortest edge is negligible. A thin but large face is not
 * one, and neither is a small equilateral face.
 */
auto
hasSliver(const MarkedPolyhedron &mesh, double maxCollapsedLength) -> bool
{
  const double maxCollapsedSquared = maxCollapsedLength * maxCollapsedLength;
  const auto   vertexPointMap      = get(CGAL::vertex_point, mesh);

  for (const auto &face : faces(mesh)) {
    double shortest = std::numeric_limits<double>::max();
    double longest  = 0.0;

    for (const auto &halfedge :
         CGAL::halfedges_around_face(CGAL::halfedge(face, mesh), mesh)) {
      const double squaredLength = CGAL::to_double(CGAL::squared_distance(
          get(vertexPointMap, CGAL::source(halfedge, mesh)),
          get(vertexPointMap, CGAL::target(halfedge, mesh))));
      shortest                   = std::min(shortest, squaredLength);
      longest                    = std::max(longest, squaredLength);
    }

    if (shortest < maxCollapsedSquared &&
        longest > NEEDLE_RATIO * NEEDLE_RATIO * shortest) {
      return true;
    }
  }
  return false;
}

} // namespace

void
collapseSlivers(MarkedPolyhedron &mesh, double maxEdgeLength)
{
  const double maxCollapsedLength =
      maxEdgeLength > 0.0 ? maxEdgeLength
                          : SLIVER_RELATIVE_LENGTH * boundingBoxDiagonal(mesh);

  // Leave a clean mesh strictly alone: no copy, no repair, no warning.
  if (!hasSliver(mesh, maxCollapsedLength)) {
    return;
  }

  MarkedPolyhedron repaired = mesh;

  // The same bound caps the flips: a face is only flattened when the triangle
  // it would replace is itself negligible.
  PMP::remove_almost_degenerate_faces(
      repaired, CGAL::parameters::needle_threshold(NEEDLE_RATIO)
                    .collapse_length_threshold(maxCollapsedLength)
                    .flip_triangle_height_threshold(maxCollapsedLength));

  // Collapsing an edge moves a vertex, which can push it through a
  // neighbouring face. Keeping the slivers is then the lesser evil: they only
  // bite downstream, whereas an intersecting shell is invalid right away.
  if (PMP::does_self_intersect(repaired)) {
    SFCGAL_WARNING("collapseSlivers: collapsing the sliver faces would make "
                   "the mesh self intersect, they are left in place");
    return;
  }

  mesh = std::move(repaired);

  if (hasSliver(mesh, maxCollapsedLength)) {
    std::ostringstream message;
    message << "collapseSlivers: some sliver faces resisted the collapse, the "
               "mesh keeps edges shorter than "
            << maxCollapsedLength
            << ". Writing it with a coarser precision will degenerate them.";
    SFCGAL_WARNING(message.str());
  }
}

} // namespace SFCGAL::detail::algorithm

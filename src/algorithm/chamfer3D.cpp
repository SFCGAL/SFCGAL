// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/chamfer3D.h"

#include "SFCGAL/algorithm/difference.h"
#include "SFCGAL/algorithm/orientedPrism.h"

namespace SFCGAL::algorithm {

auto
chamfer3D(const Geometry &solid_geom, const Geometry &edges, double distance1,
          double distance2) -> std::unique_ptr<Geometry>
{
  // 1. Generate the continuous triangular prisms (subtraction volume) along the
  // edges. The orientedPrism algorithm handles both LineString and
  // MultiLineString and manages sweep and corner connections.
  auto subtraction_geom =
      orientedPrism(solid_geom, edges, distance1, distance2);

  if (!subtraction_geom || subtraction_geom->isEmpty()) {
    return solid_geom.clone();
  }

  // 2. Perform 3D Boolean Difference to remove the prisms from the solid
  //    result = solid - (prisms)
  //
  // Using difference3D_nef which operates entirely with Nef_polyhedron_3.
  // This provides robust boolean operations with exact arithmetic and avoids
  // vertex coherence issues that occur with the MarkedPolyhedron/corefine
  // approach.
  return difference3D_nef(solid_geom, *subtraction_geom);
}

} // namespace SFCGAL::algorithm

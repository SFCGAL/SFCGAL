// Copyright (c) 2026-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/detail/algorithm/mergeCoplanarFaces.h"

#include <CGAL/Polygon_mesh_processing/compute_normal.h>

namespace PMP = CGAL::Polygon_mesh_processing;

namespace SFCGAL::algorithm::detail {

// ----------------------------------------------------------------------------------
// -- private interface
// ----------------------------------------------------------------------------------
/// @{
/// @privatesection

auto
areCoplanarFaces(const Surface_mesh_3 &mesh, faceDescriptor face1,
                 faceDescriptor face2, const Kernel::FT &epsAngle,
                 const Kernel::FT &epsDist) -> bool
{
  const Vector_3 normal1 = PMP::compute_face_normal(face1, mesh);
  const Vector_3 normal2 = PMP::compute_face_normal(face2, mesh);

  // Test if the normals are parallels
  const Kernel::FT deg2rad = CGAL_PI / Kernel::FT(180.0);
  const Kernel::FT cosEps  = std::cos(CGAL::to_double(epsAngle * deg2rad));
  const Kernel::FT dot =
      std::clamp(normal1 * normal2, Kernel::FT(-1.0), Kernel::FT(1.0));
  if (CGAL::abs(dot) < cosEps) {
    return false;
  }

  // Test if a vertex of face2 lies in the plane of face1
  const Point_3 &pt0 = mesh.point(mesh.source(mesh.halfedge(face1)));
  const CGAL::Plane_3<Kernel> plane(pt0, normal1);
  for (auto halfedge2 :
       CGAL::halfedges_around_face(mesh.halfedge(face2), mesh)) {
    const Point_3   &pt   = mesh.point(mesh.source(halfedge2));
    const Kernel::FT dist = CGAL::squared_distance(plane, pt);
    if (dist > epsDist * epsDist) {
      return false;
    }
  }

  return true;
}

auto
canJoinFace(const Surface_mesh_3 &mesh, halfedgeDescriptor halfedge) -> bool
{
  if (mesh.is_border(mesh.edge(halfedge))) {
    return false;
  }

  const faceDescriptor face         = mesh.face(halfedge);
  const faceDescriptor oppositeFace = mesh.face(mesh.opposite(halfedge));

  if (face == Surface_mesh_3::null_face() ||
      oppositeFace == Surface_mesh_3::null_face() || face == oppositeFace) {
    return false;
  }

  // If a vertex has degree < 3 before the operation,
  // removing the edge would reduce its degree below 2,
  // creating an invalid vertex.
  if (mesh.degree(mesh.source(halfedge)) < 3 ||
      mesh.degree(mesh.target(halfedge)) < 3) {
    return false;
  }

  return true;
}

/// @} end of private section

} // namespace SFCGAL::algorithm::detail

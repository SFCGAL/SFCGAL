// Copyright (c) 2026-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/mergeCoplanarFaces.h"
#include "SFCGAL/Geometry.h"
#include "SFCGAL/Kernel.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/TriangulatedSurface.h"
#include "SFCGAL/detail/algorithm/mergeCoplanarFaces.h"
#include <memory>

namespace SFCGAL::algorithm {

// ----------------------------------------------------------------------------------
// -- private interface
// ----------------------------------------------------------------------------------
/// @{
/// @privatesection

using halfedgeDescriptor = Surface_mesh_3::Halfedge_index;
using edgeDescriptor     = Surface_mesh_3::Edge_index;
using faceDescriptor     = Surface_mesh_3::Face_index;

/// @} end of private section

void
mergeCoplanarFaces(Surface_mesh_3 &mesh, const Kernel::FT &epsAngle,
                   const Kernel::FT &epsDist)
{
  std::vector<edgeDescriptor> worklist;
  worklist.reserve(mesh.num_edges());

  for (edgeDescriptor edge : mesh.edges()) {
    if (!mesh.is_border(edge)) {
      worklist.push_back(edge);
    }
  }

  while (!worklist.empty()) {
    const edgeDescriptor edge = worklist.back();
    worklist.pop_back();

    // Edge has already been removed
    if (mesh.is_removed(edge)) {
      continue;
    }

    const halfedgeDescriptor halfedge = mesh.halfedge(edge, 0);

    if (!detail::canJoinFace(mesh, halfedge)) {
      continue;
    }

    const faceDescriptor face         = mesh.face(halfedge);
    const faceDescriptor oppositeFace = mesh.face(mesh.opposite(halfedge));

    if (!detail::areCoplanarFaces(mesh, face, oppositeFace, epsAngle,
                                  epsDist)) {
      continue;
    }

    // Collect neighboring edges before the merge
    // descriptors become invalid after join_face
    std::vector<edgeDescriptor> neighbors;
    neighbors.reserve(16);
    halfedgeDescriptor neighborHalfedge = halfedge;
    do {
      neighbors.push_back(mesh.edge(neighborHalfedge));
      neighborHalfedge = mesh.next(neighborHalfedge);
    } while (neighborHalfedge != halfedge);

    neighborHalfedge = mesh.opposite(halfedge);
    do {
      neighbors.push_back(mesh.edge(neighborHalfedge));
      neighborHalfedge = mesh.next(neighborHalfedge);
    } while (neighborHalfedge != mesh.opposite(halfedge));

    // faces can be merged
    CGAL::Euler::join_face(halfedge, mesh);

    // Inspect neighbors
    for (edgeDescriptor neighborEdge : neighbors) {
      if (!mesh.is_removed(neighborEdge) && !mesh.is_border(neighborEdge)) {
        worklist.push_back(neighborEdge);
      }
    }
  }

  mesh.collect_garbage();
}

auto
mergeCoplanarFaces(const Geometry &geometry, const Kernel::FT &epsAngle,
                   const Kernel::FT &epsDist) -> std::unique_ptr<Geometry>
{
  if (geometry.isEmpty()) {
    return geometry.clone();
  }

  switch (geometry.geometryTypeId()) {
  case TYPE_POLYHEDRALSURFACE: {
    Surface_mesh_3 mesh = geometry.as<PolyhedralSurface>().toSurfaceMesh();
    mergeCoplanarFaces(mesh, epsAngle, epsDist);
    return std::make_unique<PolyhedralSurface>(mesh);
  }
  case TYPE_TRIANGULATEDSURFACE: {
    Surface_mesh_3 mesh = geometry.as<TriangulatedSurface>().toSurfaceMesh();
    mergeCoplanarFaces(mesh, epsAngle, epsDist);
    return std::make_unique<PolyhedralSurface>(mesh);
  }
  case TYPE_SOLID: {
    Surface_mesh_3 mesh = geometry.as<Solid>().toSurfaceMesh();
    mergeCoplanarFaces(mesh, epsAngle, epsDist);
    return std::make_unique<Solid>(mesh);
  }

  case TYPE_GEOMETRYCOLLECTION:
  case TYPE_LINESTRING:
  case TYPE_MULTILINESTRING:
  case TYPE_MULTIPOINT:
  case TYPE_MULTIPOLYGON:
  case TYPE_MULTISOLID:
  case TYPE_NURBSCURVE:
  case TYPE_POINT:
  case TYPE_POLYGON:
  case TYPE_TRIANGLE:
    break;
  }

  BOOST_THROW_EXCEPTION(
      Exception((boost::format("Unexpected geometry type (%s) in "
                               "SFCGAL::algorithm::mergeCoplanarFaces") %
                 geometry.geometryType())
                    .str()));
}

} // namespace SFCGAL::algorithm

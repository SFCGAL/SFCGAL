// Copyright (c) 2026-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/detail/algorithm/split3D.h"

#include "SFCGAL/Envelope.h"
#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/TriangulatedSurface.h"
#include "SFCGAL/algorithm/mergeCoplanarFaces.h"

#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/intersections.h>
#include <memory>

namespace PMP = CGAL::Polygon_mesh_processing;

namespace SFCGAL::algorithm::detail {

// ----------------------------------------------------------------------------------
// -- private interface
// ----------------------------------------------------------------------------------
/// @{
/// @privatesection

template <typename geomType>
auto
isPlaneIntersectingGeometryBBox(const CGAL::Plane_3<Kernel> &plane,
                                const geomType              &geometry) -> bool
{
  const Envelope geomEnvelope = geometry.envelope();

  // p-vertex: the corner most in the direction of the normal
  // n-vertex: the corner most opposite to the normal
  Kernel::Point_3 normalVertex(
      plane.a() >= 0 ? geomEnvelope.xMax() : geomEnvelope.xMin(),
      plane.b() >= 0 ? geomEnvelope.yMax() : geomEnvelope.yMin(),
      plane.c() >= 0 ? geomEnvelope.zMax() : geomEnvelope.zMin());
  Kernel::Point_3 oppositeVertex(
      plane.a() >= 0 ? geomEnvelope.xMin() : geomEnvelope.xMax(),
      plane.b() >= 0 ? geomEnvelope.yMin() : geomEnvelope.yMax(),
      plane.c() >= 0 ? geomEnvelope.zMin() : geomEnvelope.zMax());

  return plane.oriented_side(oppositeVertex) != CGAL::ON_POSITIVE_SIDE &&
         plane.oriented_side(normalVertex) != CGAL::ON_NEGATIVE_SIDE;
}

template <typename geomType>
auto
split3D(const geomType &geometry, const CGAL::Plane_3<Kernel> &plane,
        bool closeGeometries) -> std::unique_ptr<GeometryCollection>
{
  if (geometry.isEmpty()) {
    return std::make_unique<GeometryCollection>();
  }

  auto result = std::make_unique<GeometryCollection>();

  // if the plane does not cut the geometry bbox
  // return a collection which contains the input geometry
  if (!isPlaneIntersectingGeometryBBox(plane, geometry)) {
    result->addGeometry(geometry);
    return result;
  }

  auto workMesh = std::make_unique<Surface_mesh_3>(geometry.toSurfaceMesh());
  std::vector<std::unique_ptr<Surface_mesh_3>> meshResults;

  if (closeGeometries) {
    auto workMesh2 = std::make_unique<Surface_mesh_3>(*workMesh);
    PMP::clip(*workMesh, plane, PMP::parameters::clip_volume(true));
    PMP::clip(*workMesh2, plane.opposite(), PMP::parameters::clip_volume(true));

    // if one the mesh is empty, there is no intersection
    if (workMesh->is_empty() || workMesh2->is_empty()) {
      result->addGeometry(geometry);
      return result;
    }
    meshResults.push_back(std::move(workMesh));
    meshResults.push_back(std::move(workMesh2));
  }

  else {
    PMP::split(*workMesh, plane);
    meshResults.push_back(std::move(workMesh));
  }

  for (const auto &resultMesh : meshResults) {
    std::vector<Mesh> parts;
    PMP::split_connected_components(*resultMesh, parts);
    for (auto &splitMesh : parts) {
      mergeCoplanarFaces(splitMesh);
      if constexpr (std::is_same_v<geomType, Solid>) {
        result->addGeometry(Solid(splitMesh));
      } else {
        result->addGeometry(PolyhedralSurface(splitMesh));
      }
    }
  }

  return result;
}

// explicit instanciation
template auto
split3D<PolyhedralSurface>(const PolyhedralSurface &,
                           const CGAL::Plane_3<Kernel> &, bool)
    -> std::unique_ptr<GeometryCollection>;

template auto
split3D<TriangulatedSurface>(const TriangulatedSurface &,
                             const CGAL::Plane_3<Kernel> &, bool)
    -> std::unique_ptr<GeometryCollection>;

template auto
split3D<Solid>(const Solid &, const CGAL::Plane_3<Kernel> &, bool)
    -> std::unique_ptr<GeometryCollection>;

/// @} end of private section

} // namespace SFCGAL::algorithm::detail

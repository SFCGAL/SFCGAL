// Copyright (c) 2026-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/detail/algorithm/split.h"

#include "SFCGAL/Envelope.h"
#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/MultiSolid.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/TriangulatedSurface.h"

#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/intersections.h>

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
split(const geomType &geometry, const CGAL::Plane_3<Kernel> &plane,
      bool closeGeometries) -> std::unique_ptr<GeometryCollection>
{
  if (geometry.isEmpty()) {
    return std::make_unique<GeometryCollection>();
  }

  // if the plane does not cut the geometry bbox
  // return an empty collection
  if (!isPlaneIntersectingGeometryBBox(plane, geometry)) {
    return std::make_unique<GeometryCollection>();
  }

  auto workMesh = std::make_unique<Surface_mesh_3>(geometry.toSurfaceMesh());
  std::vector<std::unique_ptr<Surface_mesh_3>> meshResults;

  if (closeGeometries) {
    auto workMesh2 = std::make_unique<Surface_mesh_3>(*workMesh);
    PMP::clip(*workMesh, plane, PMP::parameters::clip_volume(true));
    PMP::clip(*workMesh2, plane.opposite(), PMP::parameters::clip_volume(true));

    // if one the mesh is empty, there is no intersection
    if (workMesh->is_empty() || workMesh2->is_empty()) {
      return std::make_unique<GeometryCollection>();
    }
    meshResults.push_back(std::move(workMesh));
    meshResults.push_back(std::move(workMesh2));
  }

  else {
    PMP::split(*workMesh, plane);
    meshResults.push_back(std::move(workMesh));
  }

  auto result = std::make_unique<GeometryCollection>();

  for (const auto &resultMesh : meshResults) {
    std::vector<Mesh> parts;
    PMP::split_connected_components(*resultMesh, parts);
    for (const auto &splitMesh : parts) {
      if constexpr (std::is_same_v<geomType, Solid>) {
        result->addGeometry(Solid(splitMesh));
      } else {
        result->addGeometry(PolyhedralSurface(splitMesh));
      }
    }
  }

  // if the result contains only one geometry, there is no intersection
  // return an empty result
  if (result->numGeometries() < 2) {
    return std::make_unique<GeometryCollection>();
  }

  return result;
}

template <>
auto
split(const MultiPolygon &geometry, const CGAL::Plane_3<Kernel> &plane,
      bool closeGeometries) -> std::unique_ptr<GeometryCollection>
{
  auto result = std::make_unique<GeometryCollection>();

  if (geometry.isEmpty()) {
    return result;
  }

  for (const auto &polygon : geometry) {
    std::unique_ptr<GeometryCollection> newGeoms =
        split(polygon.as<Polygon>(), plane, closeGeometries);
    for (const auto &newGeom : *newGeoms) {
      result->addGeometry(newGeom);
    }
  }

  return result;
}

template <>
auto
split(const MultiSolid &geometry, const CGAL::Plane_3<Kernel> &plane,
      bool /*unused*/) -> std::unique_ptr<GeometryCollection>
{
  auto result = std::make_unique<GeometryCollection>();

  if (geometry.isEmpty()) {
    return result;
  }

  for (const auto &solid : geometry) {
    // A solid is always closed
    std::unique_ptr<GeometryCollection> newGeoms =
        split(solid.as<Solid>(), plane, true);
    for (const auto &newGeom : *newGeoms) {
      result->addGeometry(newGeom);
    }
  }

  return result;
}

// explicit instanciation
template auto
split<PolyhedralSurface>(const PolyhedralSurface &,
                         const CGAL::Plane_3<Kernel> &, bool)
    -> std::unique_ptr<GeometryCollection>;

template auto
split<TriangulatedSurface>(const TriangulatedSurface &,
                           const CGAL::Plane_3<Kernel> &, bool)
    -> std::unique_ptr<GeometryCollection>;

/// @} end of private section

} // namespace SFCGAL::algorithm::detail

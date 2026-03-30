// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/extrude.h"

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiLineString.h"
#include "SFCGAL/MultiPoint.h"
#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/MultiSolid.h"
#include "SFCGAL/NURBSCurve.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/TriangulatedSurface.h"

#include "SFCGAL/Exception.h"

#include "SFCGAL/algorithm/force3D.h"
#include "SFCGAL/algorithm/isValid.h"
#include "SFCGAL/algorithm/normal.h"
#include "SFCGAL/algorithm/translate.h"

#include "SFCGAL/detail/tools/Log.h"

#include <utility>

namespace SFCGAL::algorithm {

// ----------------------------------------------------------------------------------
// -- private interface
// ----------------------------------------------------------------------------------
/// @{
/// @privatesection

auto
extrude(const Point &g, const Kernel::Vector_3 &extrusionVector)
    -> std::unique_ptr<LineString>;
auto
extrude(const LineString &g, const Kernel::Vector_3 &extrusionVector)
    -> PolyhedralSurface *;
auto
extrude(const Polygon &g, const Kernel::Vector_3 &extrusionVector,
        bool addTop = true) -> std::unique_ptr<Solid>;
auto
extrude(const Triangle &g, const Kernel::Vector_3 &extrusionVector)
    -> std::unique_ptr<Solid>;

auto
extrude(const MultiPoint &g, const Kernel::Vector_3 &extrusionVector)
    -> MultiLineString *;
auto
extrude(const MultiLineString &g, const Kernel::Vector_3 &extrusionVector)
    -> PolyhedralSurface *;
auto
extrude(const MultiPolygon &g, const Kernel::Vector_3 &extrusionVector)
    -> MultiSolid *;

auto
extrude(const TriangulatedSurface &g, const Kernel::Vector_3 &extrusionVector)
    -> Solid *;
auto
extrude(const PolyhedralSurface &g, const Kernel::Vector_3 &extrusionVector)
    -> Solid *;

auto
extrude(const GeometryCollection &g, const Kernel::Vector_3 &extrusionVector)
    -> GeometryCollection *;

auto
extrude(const Point &g, const Kernel::Vector_3 &extrusionVector)
    -> std::unique_ptr<LineString>
{
  if (g.isEmpty()) {
    return std::make_unique<LineString>();
  }

  Kernel::Point_3 const a = g.toPoint_3();
  Kernel::Point_3 const b = a + extrusionVector;

  return std::make_unique<LineString>(Point(a), Point(b));
}

auto
extrude(const LineString &g, const Kernel::Vector_3 &extrusionVector)
    -> PolyhedralSurface *
{

  std::unique_ptr<PolyhedralSurface> polyhedralSurface(new PolyhedralSurface());

  if (g.isEmpty()) {
    return polyhedralSurface.release();
  }

  for (size_t i = 0; i < g.numPoints() - 1; i++) {
    std::unique_ptr<LineString> ring(new LineString);

    Kernel::Point_3 const a = g.pointN(i).toPoint_3();
    Kernel::Point_3 const b = g.pointN(i + 1).toPoint_3();
    ring->addPoint(std::make_unique<Point>(a));
    ring->addPoint(std::make_unique<Point>(b));
    ring->addPoint(std::make_unique<Point>(b + extrusionVector));
    ring->addPoint(std::make_unique<Point>(a + extrusionVector));
    ring->addPoint(std::make_unique<Point>(a));

    polyhedralSurface->addPatch(std::make_unique<Polygon>(std::move(ring)));
  }

  return polyhedralSurface.release();
}

auto
extrude(const Polygon &g, const Kernel::Vector_3 &extrusionVector, bool addTop)
    -> std::unique_ptr<Solid>
{
  if (g.isEmpty()) {
    return std::make_unique<Solid>();
  }

  bool const reverseOrientation = (extrusionVector * normal3D<Kernel>(g)) > 0;

  // resulting shell
  PolyhedralSurface polyhedralSurface;

  // "bottom"
  Polygon bottom(g);
  force3D(bottom);

  if (reverseOrientation) {
    bottom.reverse();
  }

  polyhedralSurface.addPatch(bottom);

  // "top"
  if (addTop) {
    Polygon top(bottom);
    top.reverse();
    translate(top, extrusionVector);
    polyhedralSurface.addPatch(top);
  }
  // exterior ring and interior rings extruded
  for (size_t i = 0; i < bottom.numRings(); i++) {
    std::unique_ptr<PolyhedralSurface> boundaryExtruded(
        extrude(bottom.ringN(i), extrusionVector));

    for (size_t j = 0; j < boundaryExtruded->numPatches(); j++) {
      boundaryExtruded->patchN(j).reverse();
      polyhedralSurface.addPatch(boundaryExtruded->patchN(j));
    }
  }

  return std::make_unique<Solid>(polyhedralSurface);
}

auto
extrude(const Triangle &g, const Kernel::Vector_3 &extrusionVector)
    -> std::unique_ptr<Solid>
{
  return extrude(g.toPolygon(), extrusionVector);
}

auto
extrude(const MultiPoint &g, const Kernel::Vector_3 &extrusionVector)
    -> MultiLineString *
{
  std::unique_ptr<MultiLineString> result(new MultiLineString());

  if (g.isEmpty()) {
    return result.release();
  }

  for (size_t i = 0; i < g.numGeometries(); i++) {
    result->addGeometry(extrude(g.pointN(i), extrusionVector));
  }

  return result.release();
}

auto
extrude(const MultiLineString &g, const Kernel::Vector_3 &extrusionVector)
    -> PolyhedralSurface *
{
  std::unique_ptr<PolyhedralSurface> result(new PolyhedralSurface());

  if (g.isEmpty()) {
    return result.release();
  }

  for (size_t i = 0; i < g.numGeometries(); i++) {
    std::unique_ptr<PolyhedralSurface> extruded(
        extrude(g.lineStringN(i), extrusionVector));

    for (size_t j = 0; j < extruded->numPatches(); j++) {
      result->addPatch(extruded->patchN(j));
    }
  }

  return result.release();
}

auto
extrude(const MultiPolygon &g, const Kernel::Vector_3 &extrusionVector)
    -> MultiSolid *
{
  std::unique_ptr<MultiSolid> result(new MultiSolid());

  if (g.isEmpty()) {
    return result.release();
  }

  for (size_t i = 0; i < g.numGeometries(); i++) {
    result->addGeometry(extrude(g.polygonN(i), extrusionVector));
  }

  return result.release();
}

auto
extrude(const TriangulatedSurface &g, const Kernel::Vector_3 &extrusionVector)
    -> Solid *
{
  std::unique_ptr<Solid> result(new Solid());

  if (g.isEmpty()) {
    return result.release();
  }

  // bottom and top
  for (size_t i = 0; i < g.numPatches(); i++) {
    Triangle bottomPart(g.patchN(i));
    force3D(bottomPart);
    bottomPart.reverse();
    result->exteriorShell().addPatch(bottomPart);

    Triangle topPart(g.patchN(i));
    force3D(topPart);
    translate(topPart, extrusionVector);
    result->exteriorShell().addPatch(topPart);
  }

  // boundary
  std::unique_ptr<Geometry> boundary(g.boundary());
  BOOST_ASSERT(boundary.get() != NULL);

  // closed surface extruded
  if (!boundary->isEmpty()) {
    std::unique_ptr<Geometry> extrudedBoundary(
        extrude(*boundary, extrusionVector));
    BOOST_ASSERT(extrudedBoundary->is<PolyhedralSurface>());
    result->exteriorShell().addPatches(
        extrudedBoundary->as<PolyhedralSurface>());
  }

  return result.release();
}

auto
extrude(const PolyhedralSurface &g, const Kernel::Vector_3 &extrusionVector)
    -> Solid *
{
  if (g.isEmpty()) {
    return new Solid();
  }

  TriangulatedSurface triangulatedSurface;
  triangulate::triangulatePolygon3D(g, triangulatedSurface);
  return extrude(triangulatedSurface, extrusionVector);
}

auto
extrude(const GeometryCollection &g, const Kernel::Vector_3 &extrusionVector)
    -> GeometryCollection *
{
  std::unique_ptr<GeometryCollection> result(new GeometryCollection());

  if (g.isEmpty()) {
    return result.release();
  }

  for (size_t i = 0; i < g.numGeometries(); i++) {
    result->addGeometry(extrude(g.geometryN(i), extrusionVector));
  }

  return result.release();
}

/// @private
auto
extrude(const Geometry &inputGeometry, const Kernel::Vector_3 &extrusionVector)
    -> std::unique_ptr<Geometry>
{
  switch (inputGeometry.geometryTypeId()) {
  case TYPE_POINT:
    return std::unique_ptr<Geometry>(
        extrude(inputGeometry.as<Point>(), extrusionVector));

  case TYPE_LINESTRING:
    return std::unique_ptr<Geometry>(
        extrude(inputGeometry.as<LineString>(), extrusionVector));

  case TYPE_NURBSCURVE: {
    auto lineString =
        inputGeometry.as<NURBSCurve>().toLineString(); // default parameters
    if (!lineString || lineString->isEmpty()) {
      return std::make_unique<PolyhedralSurface>(); // empty result
    }
    return std::unique_ptr<Geometry>(extrude(*lineString, extrusionVector));
  }

  case TYPE_POLYGON:
    return std::unique_ptr<Geometry>(
        extrude(inputGeometry.as<Polygon>(), extrusionVector));

  case TYPE_TRIANGLE:
    return std::unique_ptr<Geometry>(
        extrude(inputGeometry.as<Triangle>(), extrusionVector));

  case TYPE_GEOMETRYCOLLECTION:
    return std::unique_ptr<Geometry>(
        extrude(inputGeometry.as<GeometryCollection>(), extrusionVector));

  case TYPE_MULTIPOINT:
    return std::unique_ptr<Geometry>(
        extrude(inputGeometry.as<MultiPoint>(), extrusionVector));

  case TYPE_MULTILINESTRING:
    return std::unique_ptr<Geometry>(
        extrude(inputGeometry.as<MultiLineString>(), extrusionVector));

  case TYPE_MULTIPOLYGON:
    return std::unique_ptr<Geometry>(
        extrude(inputGeometry.as<MultiPolygon>(), extrusionVector));

  case TYPE_TRIANGULATEDSURFACE:
    return std::unique_ptr<Geometry>(
        extrude(inputGeometry.as<TriangulatedSurface>(), extrusionVector));

  case TYPE_POLYHEDRALSURFACE:
    return std::unique_ptr<Geometry>(
        extrude(inputGeometry.as<PolyhedralSurface>(), extrusionVector));

  case TYPE_SOLID:
  case TYPE_MULTISOLID:
    // extrusion not available
    break;
  }

  BOOST_THROW_EXCEPTION(InappropriateGeometryException(
      (boost::format("Extrusion of %s is not supported") %
       inputGeometry.geometryType())
          .str()));
}

/// @private
auto
extrude(const Geometry &inputGeom, const Kernel::FT &displacementX,
        const Kernel::FT &displacementY, const Kernel::FT &displacementZ,
        NoValidityCheck /*unused*/) -> std::unique_ptr<Geometry>
{
  return extrude(inputGeom,
                 Kernel::Vector_3(displacementX, displacementY, displacementZ));
}

/// @private
auto
extrude(const Geometry &geometry, const Kernel::FT &deltaX,
        const Kernel::FT &deltaY, const Kernel::FT &deltaZ)
    -> std::unique_ptr<Geometry>
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY(geometry);
  std::unique_ptr<Geometry> result(
      extrude(geometry, deltaX, deltaY, deltaZ, NoValidityCheck()));
  propagateValidityFlag(*result, true);
  return result;
}

/// @private
SFCGAL_API auto
extrude(const Geometry &geom, const double &displacementX,
        const double &displacementY, const double &displacementZ)
    -> std::unique_ptr<Geometry>
{
  if (!std::isfinite(displacementX) || !std::isfinite(displacementY) ||
      !std::isfinite(displacementZ)) {
    BOOST_THROW_EXCEPTION(NonFiniteValueException(
        "trying to extrude with non finite value in extrusionVector"));
  }

  return extrude(geom, Kernel::FT(displacementX), Kernel::FT(displacementY),
                 Kernel::FT(displacementZ));
}

/// @private
SFCGAL_API auto
extrude(const Polygon &polygon, const double &height)
    -> std::unique_ptr<Geometry>
{

  if (!std::isfinite(height)) {
    BOOST_THROW_EXCEPTION(NonFiniteValueException(
        "trying to extrude with non finite value in extrusionVector"));
  }

  return std::unique_ptr<Geometry>(
      extrude(polygon, Kernel::Vector_3(0.0, 0.0, height), false));
}
} // namespace SFCGAL::algorithm

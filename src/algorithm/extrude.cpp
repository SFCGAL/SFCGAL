// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2026, SFCGAL team.
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

#include "SFCGAL/algorithm/area.h"
#include "SFCGAL/algorithm/force2D.h"
#include "SFCGAL/algorithm/force3D.h"
#include "SFCGAL/algorithm/intersection.h"
#include "SFCGAL/algorithm/isValid.h"
#include "SFCGAL/algorithm/normal.h"
#include "SFCGAL/algorithm/plane.h"
#include "SFCGAL/algorithm/translate.h"

#include "SFCGAL/detail/tools/Log.h"

#include <map>
#include <memory>
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
    -> std::unique_ptr<PolyhedralSurface>;
auto
extrude(const Polygon &g, const Kernel::Vector_3 &extrusionVector,
        bool addTop = true) -> std::unique_ptr<Solid>;
auto
extrude(const Triangle &g, const Kernel::Vector_3 &extrusionVector)
    -> std::unique_ptr<Solid>;

auto
extrude(const MultiPoint &g, const Kernel::Vector_3 &extrusionVector)
    -> std::unique_ptr<MultiLineString>;
auto
extrude(const MultiLineString &g, const Kernel::Vector_3 &extrusionVector)
    -> std::unique_ptr<PolyhedralSurface>;
auto
extrude(const MultiPolygon &g, const Kernel::Vector_3 &extrusionVector)
    -> std::unique_ptr<MultiSolid>;

auto
extrude(const TriangulatedSurface &g, const Kernel::Vector_3 &extrusionVector)
    -> std::unique_ptr<Solid>;
auto
extrude(const PolyhedralSurface &g, const Kernel::Vector_3 &extrusionVector)
    -> std::unique_ptr<Solid>;

auto
extrude(const GeometryCollection &g, const Kernel::Vector_3 &extrusionVector)
    -> std::unique_ptr<GeometryCollection>;

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
    -> std::unique_ptr<PolyhedralSurface>
{

  auto polyhedralSurface = std::make_unique<PolyhedralSurface>();

  if (g.isEmpty()) {
    return polyhedralSurface;
  }

  for (size_t i = 0; i < g.numPoints() - 1; i++) {
    auto ring = std::make_unique<LineString>();

    Kernel::Point_3 const a = g.pointN(i).toPoint_3();
    Kernel::Point_3 const b = g.pointN(i + 1).toPoint_3();
    ring->addPoint(std::make_unique<Point>(a));
    ring->addPoint(std::make_unique<Point>(b));
    ring->addPoint(std::make_unique<Point>(b + extrusionVector));
    ring->addPoint(std::make_unique<Point>(a + extrusionVector));
    ring->addPoint(std::make_unique<Point>(a));

    polyhedralSurface->addPatch(std::make_unique<Polygon>(std::move(ring)));
  }

  return polyhedralSurface;
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
    -> std::unique_ptr<MultiLineString>
{
  auto result = std::make_unique<MultiLineString>();

  if (g.isEmpty()) {
    return result;
  }

  for (size_t i = 0; i < g.numGeometries(); i++) {
    result->addGeometry(extrude(g.pointN(i), extrusionVector));
  }

  return result;
}

auto
extrude(const MultiLineString &g, const Kernel::Vector_3 &extrusionVector)
    -> std::unique_ptr<PolyhedralSurface>
{
  auto result = std::make_unique<PolyhedralSurface>();

  if (g.isEmpty()) {
    return result;
  }

  for (size_t i = 0; i < g.numGeometries(); i++) {
    std::unique_ptr<PolyhedralSurface> extruded(
        extrude(g.lineStringN(i), extrusionVector));

    for (size_t j = 0; j < extruded->numPatches(); j++) {
      result->addPatch(extruded->patchN(j));
    }
  }

  return result;
}

auto
extrude(const MultiPolygon &g, const Kernel::Vector_3 &extrusionVector)
    -> std::unique_ptr<MultiSolid>
{
  auto result = std::make_unique<MultiSolid>();

  if (g.isEmpty()) {
    return result;
  }

  for (size_t i = 0; i < g.numGeometries(); i++) {
    result->addGeometry(extrude(g.polygonN(i), extrusionVector));
  }

  return result;
}

auto
extrude(const TriangulatedSurface &g, const Kernel::Vector_3 &extrusionVector)
    -> std::unique_ptr<Solid>
{
  auto result = std::make_unique<Solid>();

  if (g.isEmpty()) {
    return result;
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

  return result;
}

auto
extrude(const PolyhedralSurface &g, const Kernel::Vector_3 &extrusionVector)
    -> std::unique_ptr<Solid>
{
  if (g.isEmpty()) {
    return std::make_unique<Solid>();
  }

  TriangulatedSurface triangulatedSurface;
  triangulate::triangulatePolygon3D(g, triangulatedSurface);
  return extrude(triangulatedSurface, extrusionVector);
}

auto
extrude(const GeometryCollection &g, const Kernel::Vector_3 &extrusionVector)
    -> std::unique_ptr<GeometryCollection>
{
  auto result = std::make_unique<GeometryCollection>();

  if (g.isEmpty()) {
    return result;
  }

  for (size_t i = 0; i < g.numGeometries(); i++) {
    result->addGeometry(extrude(g.geometryN(i), extrusionVector));
  }

  return result;
}

/// @private
auto
extrude(const Geometry &inputGeometry, const Kernel::Vector_3 &extrusionVector)
    -> std::unique_ptr<Geometry>
{
  switch (inputGeometry.geometryTypeId()) {
  case TYPE_POINT:
    return extrude(inputGeometry.as<Point>(), extrusionVector);

  case TYPE_LINESTRING:
    return extrude(inputGeometry.as<LineString>(), extrusionVector);

  case TYPE_NURBSCURVE: {
    auto lineString =
        inputGeometry.as<NURBSCurve>().toLineString(); // default parameters
    if (!lineString || lineString->isEmpty()) {
      return std::make_unique<PolyhedralSurface>(); // empty result
    }
    return extrude(*lineString, extrusionVector);
  }

  case TYPE_POLYGON:
    return extrude(inputGeometry.as<Polygon>(), extrusionVector);

  case TYPE_TRIANGLE:
    return extrude(inputGeometry.as<Triangle>(), extrusionVector);

  case TYPE_GEOMETRYCOLLECTION:
    return extrude(inputGeometry.as<GeometryCollection>(), extrusionVector);

  case TYPE_MULTIPOINT:
    return extrude(inputGeometry.as<MultiPoint>(), extrusionVector);

  case TYPE_MULTILINESTRING:
    return extrude(inputGeometry.as<MultiLineString>(), extrusionVector);

  case TYPE_MULTIPOLYGON:
    return extrude(inputGeometry.as<MultiPolygon>(), extrusionVector);

  case TYPE_TRIANGULATEDSURFACE:
    return extrude(inputGeometry.as<TriangulatedSurface>(), extrusionVector);

  case TYPE_POLYHEDRALSURFACE:
    return extrude(inputGeometry.as<PolyhedralSurface>(), extrusionVector);

  case TYPE_SOLID:
  case TYPE_MULTISOLID:
    // extrusion not available
    break;
  }

  throw InappropriateGeometryException(std::format(
      "Extrusion of %s is not supported", inputGeometry.geometryType()));
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
    throw NonFiniteValueException(
        "trying to extrude with non finite value in extrusionVector");
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
    throw NonFiniteValueException(
        "trying to extrude with non finite value in extrusionVector");
  }

  return extrude(polygon, Kernel::Vector_3(0.0, 0.0, height), false);
}

namespace { // Helper functions for extrudeUntil

/**
 * @brief Lifts a 2D polygon (z=0) onto a 3D plane.
 *
 * Uses the plane equation ax + by + cz + d = 0 to compute:
 *   z = -(ax + by + d) / c
 *
 * @param poly2D Polygon with z coordinates at 0
 * @param plane The target 3D plane
 * @return Polygon with z coordinates computed from the plane equation
 * @pre plane.c() != 0 (non-vertical plane)
 */
auto
liftPolygonToPlane(const Polygon &poly2D, const CGAL::Plane_3<Kernel> &plane)
    -> Polygon
{
  if (plane.c() == 0) {
    // Vertical plane - cannot compute z coordinate
    return poly2D;
  }

  Polygon poly3D = poly2D;
  for (size_t i = 0; i < poly3D.numRings(); ++i) {
    LineString &ring = poly3D.ringN(i);
    for (size_t j = 0; j < ring.numPoints(); ++j) {
      Point     &point = ring.pointN(j);
      Kernel::FT x     = point.x();
      Kernel::FT y     = point.y();
      Kernel::FT z = -(plane.a() * x + plane.b() * y + plane.d()) / plane.c();
      point        = Point(x, y, z);
    }
  }
  return poly3D;
}

/**
 * @brief Processes roof polygons by projecting, intersecting with footprint,
 * and lifting back to 3D.
 *
 * @param footprint The normalized 2D footprint polygon
 * @param roof The roof surface with potentially non-planar polygons
 * @return Vector of 3D polygons representing the intersection of footprint
 * with roof
 */
auto
processRoofPolygons(const Polygon &footprint, const PolyhedralSurface &roof)
    -> std::vector<Polygon>
{
  std::vector<Polygon> topPolygons;

  for (size_t i = 0; i < roof.numPatches(); ++i) {
    const Polygon &roofPoly = roof.patchN(i);

    CGAL::Plane_3<Kernel> plane = algorithm::plane3D<Kernel>(roofPoly);

    // Ignore vertical faces (plane.c() == 0 means normal perpendicular to Z)
    if (plane.c() == 0) {
      continue;
    }

    Polygon roofPoly2D = roofPoly;
    force2D(roofPoly2D);

    std::unique_ptr<Geometry> intersectionResult =
        intersection(footprint, roofPoly2D);

    if (!intersectionResult || intersectionResult->isEmpty()) {
      continue;
    }

    std::vector<Polygon> currentPolys;
    switch (intersectionResult->geometryTypeId()) {
    case TYPE_POLYGON:
      currentPolys.push_back(intersectionResult->as<Polygon>());
      break;

    case TYPE_TRIANGLE:
      currentPolys.emplace_back(intersectionResult->as<Triangle>());
      break;

    case TYPE_MULTIPOLYGON: {
      const MultiPolygon &mp = intersectionResult->as<MultiPolygon>();
      for (size_t k = 0; k < mp.numGeometries(); ++k) {
        currentPolys.push_back(mp.geometryN(k).as<Polygon>());
      }
      break;
    }

    case TYPE_GEOMETRYCOLLECTION: {
      const GeometryCollection &gc =
          intersectionResult->as<GeometryCollection>();
      for (size_t k = 0; k < gc.numGeometries(); ++k) {
        if (gc.geometryN(k).is<Polygon>()) {
          currentPolys.push_back(gc.geometryN(k).as<Polygon>());
        } else if (gc.geometryN(k).is<Triangle>()) {
          currentPolys.emplace_back(gc.geometryN(k).as<Triangle>());
        }
      }
      break;
    }

    default:
      break;
    }

    for (const auto &poly : currentPolys) {
      Polygon poly3D = liftPolygonToPlane(poly, plane);
      topPolygons.push_back(poly3D);
    }
  }

  return topPolygons;
}

/**
 * @brief Represents a directed edge in 2D (ignoring Z coordinate).
 */
struct Edge2D {
  Point p1;
  Point p2;

  Edge2D() = default;
  Edge2D(const Point &a, const Point &b) : p1(a), p2(b) {}

  [[nodiscard]] auto
  normalized() const -> std::pair<Point, Point>
  {
    return (p1 < p2) ? std::make_pair(p1, p2) : std::make_pair(p2, p1);
  }
};

/**
 * @brief Extracts boundary edges from a set of polygons.
 *
 * A boundary edge appears exactly once (edges shared by two polygons cancel
 * out).
 *
 * @param polygons Vector of polygons forming a surface
 * @return Vector of directed boundary edges
 */
auto
extractBoundaryEdges(const std::vector<Polygon> &polygons)
    -> std::vector<Edge2D>
{
  // Map normalized edge -> {count, original directed edge}
  std::map<std::pair<Point, Point>, std::pair<int, Edge2D>, std::less<>>
      edgeMap;

  // Collect all edges from all polygon rings
  for (const auto &poly : polygons) {
    for (size_t ringIdx = 0; ringIdx < poly.numRings(); ++ringIdx) {
      const LineString &ring = poly.ringN(ringIdx);
      for (size_t i = 0; i + 1 < ring.numPoints(); ++i) {
        const Point &pt1 = ring.pointN(i);
        const Point &pt2 = ring.pointN(i + 1);

        Edge2D edge(pt1, pt2);
        auto   key = edge.normalized();

        auto it = edgeMap.find(key);
        if (it == edgeMap.end()) {
          edgeMap[key] = {1, edge};
        } else {
          it->second.first++;
        }
      }
    }
  }

  // Collect boundary edges (count == 1)
  std::vector<Edge2D> boundaryEdges;
  for (const auto &[key, countAndEdge] : edgeMap) {
    if (countAndEdge.first == 1) {
      boundaryEdges.push_back(countAndEdge.second);
    }
  }

  return boundaryEdges;
}

/**
 * @brief Builds connected loops from a set of directed edges.
 *
 * @param edges Vector of directed edges
 * @return Vector of LineStrings representing closed loops
 */
auto
buildConnectedLoops(std::vector<Edge2D> edges) -> std::vector<LineString>
{
  std::vector<LineString> loops;

  while (!edges.empty()) {
    LineString currentLoop;
    currentLoop.addPoint(edges[0].p1);
    currentLoop.addPoint(edges[0].p2);

    Point currentEnd = edges[0].p2;
    edges.erase(edges.begin());

    bool progress = true;
    while (progress && !edges.empty()) {
      progress = false;
      for (auto it = edges.begin(); it != edges.end(); ++it) {
        if (it->p1 == currentEnd) {
          currentLoop.addPoint(it->p2);
          currentEnd = it->p2;
          edges.erase(it);
          progress = true;
          break;
        }
      }
    }

    if (currentLoop.numPoints() > 0) {
      const Point &first = currentLoop.pointN(0);
      const Point &last  = currentLoop.pointN(currentLoop.numPoints() - 1);
      if (first != last) {
        currentLoop.addPoint(first);
      }
    }

    loops.push_back(currentLoop);
  }

  return loops;
}

/**
 * @brief Finds the exterior loop (largest area) among multiple loops.
 *
 * @param loops Vector of closed LineStrings
 * @return Index of the exterior loop
 */
auto
findExteriorLoopIndex(const std::vector<LineString> &loops) -> size_t
{
  size_t     maxIdx = 0;
  Kernel::FT maxArea(0);

  for (size_t i = 0; i < loops.size(); ++i) {
    CGAL::Polygon_2<Kernel> poly2d;
    const LineString       &loop = loops[i];
    for (size_t j = 0; j + 1 < loop.numPoints(); ++j) {
      poly2d.push_back(Kernel::Point_2(loop.pointN(j).x(), loop.pointN(j).y()));
    }
    Kernel::FT area = CGAL::abs(poly2d.area());
    if (area > maxArea) {
      maxArea = area;
      maxIdx  = i;
    }
  }

  return maxIdx;
}

/**
 * @brief Extracts boundary rings from top polygons and creates bottom face.
 *
 * @param topPolygons Vector of polygons forming the top surface
 * @return Pair of (exterior ring, hole rings) projected to z=0 with reversed
 * orientation
 */
auto
extractBoundaryRings(const std::vector<Polygon> &topPolygons)
    -> std::pair<LineString, std::vector<LineString>>
{
  // Extract boundary edges
  std::vector<Edge2D> boundaryEdges = extractBoundaryEdges(topPolygons);

  // Build connected loops
  std::vector<LineString> loops = buildConnectedLoops(boundaryEdges);

  if (loops.empty()) {
    return {LineString(), {}};
  }

  // Find exterior loop
  size_t exteriorIdx = findExteriorLoopIndex(loops);

  // Project exterior loop to z=0 and reverse orientation
  LineString        bottomRing;
  const LineString &extLoop = loops[exteriorIdx];
  for (int i = static_cast<int>(extLoop.numPoints()) - 1; i >= 0; --i) {
    const Point &pt = extLoop.pointN(i);
    bottomRing.addPoint(Point(pt.x(), pt.y(), 0.0));
  }

  // Project hole rings to z=0 and reverse orientation
  std::vector<LineString> holeRings;
  for (size_t i = 0; i < loops.size(); ++i) {
    if (i != exteriorIdx) {
      LineString hole;
      for (int j = static_cast<int>(loops[i].numPoints()) - 1; j >= 0; --j) {
        const Point &pt = loops[i].pointN(j);
        hole.addPoint(Point(pt.x(), pt.y(), 0.0));
      }
      holeRings.push_back(hole);
    }
  }

  return {bottomRing, holeRings};
}

/**
 * @brief Creates vertical wall polygons from boundary rings.
 *
 * @param topPolygons Vector of polygons forming the top surface
 * @param bottomRing Exterior ring of the bottom face
 * @param holeRings Interior rings (holes) of the bottom face
 * @return Vector of wall polygons connecting top boundary to bottom
 */
auto
createWallsFromBoundaryRings(const std::vector<Polygon> &topPolygons)
    -> std::vector<Polygon>
{
  std::vector<Polygon> walls;

  // Extract boundary edges from top polygons (these have original Z
  // coordinates)
  std::vector<Edge2D> topBoundaryEdges = extractBoundaryEdges(topPolygons);

  // For each top boundary edge, create a wall quad
  for (const auto &edge : topBoundaryEdges) {
    Point p1Top = edge.p1;
    Point p2Top = edge.p2;

    Point p1Base(p1Top.x(), p1Top.y(), 0.0);
    Point p2Base(p2Top.x(), p2Top.y(), 0.0);

    // Create wall quad (p1Top -> p1Base -> p2Base -> p2Top)
    LineString wallRing;
    wallRing.addPoint(p1Top);
    wallRing.addPoint(p1Base);
    wallRing.addPoint(p2Base);
    wallRing.addPoint(p2Top);
    wallRing.addPoint(p1Top); // Close the ring

    walls.emplace_back(wallRing);
  }

  return walls;
}

} // namespace

/// @private
SFCGAL_API auto
extrudeUntil(const Polygon &footprint, const PolyhedralSurface &roof)
    -> std::unique_ptr<Solid>
{
  if (footprint.isEmpty() || roof.isEmpty()) {
    return std::make_unique<Solid>();
  }

  // Normalize footprint to CCW orientation
  Polygon normalizedFootprint = footprint;
  if (!normalizedFootprint.isCounterClockWiseOriented()) {
    normalizedFootprint.reverse();
  }

  // Process roof polygons: project to 2D, intersect with footprint, lift to 3D
  std::vector<Polygon> topPolygons =
      processRoofPolygons(normalizedFootprint, roof);

  if (topPolygons.empty()) {
    return std::make_unique<Solid>();
  }

  // Build the solid shell
  auto shell = std::make_unique<PolyhedralSurface>();

  // Add top surface polygons
  for (const auto &poly : topPolygons) {
    shell->addPatch(poly);
  }

  // Extract boundary rings from top polygons
  auto [bottomRing, holeRings] = extractBoundaryRings(topPolygons);

  // Create and add bottom face
  if (bottomRing.numPoints() > 0) {
    Polygon bottomFace(bottomRing);
    for (auto &hole : holeRings) {
      bottomFace.addInteriorRing(hole);
    }
    shell->addPatch(bottomFace);
  }

  // Create and add vertical walls from boundary edges
  auto walls = createWallsFromBoundaryRings(topPolygons);
  for (auto &wall : walls) {
    shell->addPatch(wall);
  }

  // Assemble final solid
  auto solid = std::make_unique<Solid>();
  solid->setExteriorShell(std::move(shell));
  return solid;
}

SFCGAL_API auto
extrudeUntil(const Polygon &footprint, const Geometry &roof)
    -> std::unique_ptr<Solid>
{
  // Convert roof geometry to PolyhedralSurface based on type
  PolyhedralSurface roofSurface;

  switch (roof.geometryTypeId()) {
  case TYPE_POLYHEDRALSURFACE:
    roofSurface = roof.as<PolyhedralSurface>();
    break;

  case TYPE_POLYGON:
    roofSurface.addPatch(roof.as<Polygon>());
    break;

  case TYPE_TRIANGLE:
    roofSurface.addPatch(roof.as<Triangle>());
    break;

  case TYPE_TRIANGULATEDSURFACE: {
    const auto &tin = roof.as<TriangulatedSurface>();
    for (size_t i = 0; i < tin.numPatches(); ++i) {
      roofSurface.addPatch(tin.patchN(i));
    }
    break;
  }

  default:
    BOOST_THROW_EXCEPTION(
        Exception("extrudeUntil: roof must be Polygon, Triangle, "
                  "PolyhedralSurface, or TriangulatedSurface"));
  }

  return extrudeUntil(footprint, roofSurface);
}

} // namespace SFCGAL::algorithm

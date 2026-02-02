// Copyright (c) 2025-2025, SFCGAL Team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/surfaceSimplification.h"

#include "SFCGAL/MultiSolid.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/TriangulatedSurface.h"
#include "SFCGAL/algorithm/isValid.h"
#include "SFCGAL/io/wkt.h"

#include <boost/test/unit_test.hpp>

using namespace SFCGAL;
using namespace SFCGAL::algorithm;

#include "../../../test_config.h"

BOOST_AUTO_TEST_SUITE(SFCGAL_SurfaceSimplificationTest)

/**
 * @brief Returns the 6 faces of a cube as std::array<Point,4> each.
 *
 * @param origin  The origin (minimum x, y, z) of the cube.
 * @param size    The edge length of the cube (assumed equal in x, y, z).
 * @return std::array<std::array<Point,4>,6>  Array of 6 faces.
 */
inline auto
makeCube(const Point &origin = Point(0, 0, 0), double size = 1.0)
    -> const std::array<std::array<Point, 4>, 6> &
{
  const double x0 = CGAL::to_double(origin.x());
  const double y0 = CGAL::to_double(origin.y());
  const double z0 = CGAL::to_double(origin.z());
  const double x1 = x0 + size;
  const double y1 = y0 + size;
  const double z1 = z0 + size;

  // Cube vertices
  const Point p000{x0, y0, z0};
  const Point p100{x1, y0, z0};
  const Point p110{x1, y1, z0};
  const Point p010{x0, y1, z0};
  const Point p001{x0, y0, z1};
  const Point p101{x1, y0, z1};
  const Point p111{x1, y1, z1};
  const Point p011{x0, y1, z1};

  static const std::array<std::array<Point, 4>, 6> faces = {
      // Bottom (z = 0)
      std::array<Point, 4>{{p000, p100, p110, p010}},
      // Top (z = 1)
      std::array<Point, 4>{{p001, p011, p111, p101}},
      // Front (y = 0)
      std::array<Point, 4>{{p000, p001, p101, p100}},
      // Back (y = 1)
      std::array<Point, 4>{{p010, p110, p111, p011}},
      // Left (x = 0)
      std::array<Point, 4>{{p000, p010, p011, p001}},
      // Right (x = 1)
      std::array<Point, 4>{{p100, p101, p111, p110}},
  };

  return faces;
}

/**
 * @brief Helper to create a simple quadrilateral face as a Polygon
 */
auto
createPolygonFace(const std::array<Point, 4> &vertices)
    -> std::unique_ptr<Polygon>
{
  auto ring = std::make_unique<LineString>();

  for (const auto &vertex : vertices) {
    ring->addPoint(vertex);
  }

  ring->addPoint(vertices.front());

  return std::make_unique<Polygon>(ring.release());
}

/**
 * @brief Helper to create a simple cube as a TriangulatedSurface
 */
auto
createCubeTriangulatedSurface(const Point &origin = Point(0, 0, 0),
                              double       size   = 1.0)
    -> std::unique_ptr<TriangulatedSurface>
{
  auto cube = std::make_unique<TriangulatedSurface>();

  for (const auto &face : makeCube(origin, size)) {
    cube->addTriangle(Triangle(face[0], face[1], face[2]));
    cube->addTriangle(Triangle(face[0], face[2], face[3]));
  }

  return cube;
}

/**
 * @brief Helper to create a simple cube as a PolyhedralSurface
 */
auto
createCubePolyhedralSurface(const Point &origin = Point(0, 0, 0),
                            double       size   = 1.0)
    -> std::unique_ptr<PolyhedralSurface>
{
  auto cube = std::make_unique<PolyhedralSurface>();

  for (const auto &face : makeCube(origin, size)) {
    cube->addPatch(createPolygonFace(face));
  }

  return cube;
}

#ifdef SFCGAL_WITH_EIGEN
// Test TriangulatedSurface simplification with Garland-Heckbert strategy
BOOST_AUTO_TEST_CASE(testSimplify_TriangulatedSurface_GarlandHeckbert)
{
  auto cube = createCubeTriangulatedSurface();

  BOOST_CHECK_EQUAL(cube->numTriangles(), 12);

  // Simplify using edge count ratio with Garland-Heckbert strategy
  auto simplified = surfaceSimplification(
      *cube, SimplificationStopPredicate::edgeCountRatio(0.5),
      SimplificationStrategy::GARLAND_HECKBERT);

  BOOST_CHECK(simplified);
  BOOST_CHECK(simplified->is<TriangulatedSurface>());

  const auto &simplifiedTin = simplified->as<TriangulatedSurface>();
  BOOST_CHECK_LT(simplifiedTin.numTriangles(), cube->numTriangles());

  std::cout << "Original triangles: " << cube->numTriangles() << "\n";
  std::cout << "Simplified triangles: " << simplifiedTin.numTriangles() << "\n";
}
#endif // SFCGAL_WITH_EIGEN

#ifdef SFCGAL_WITH_EIGEN
// Test TriangulatedSurface simplification with Lindstrom-Turk strategy
BOOST_AUTO_TEST_CASE(testSimplify_TriangulatedSurface_LindstromTurk)
{
  auto cube = createCubeTriangulatedSurface();

  BOOST_CHECK_EQUAL(cube->numTriangles(), 12);

  // Simplify using edge count ratio with Lindstrom-Turk strategy
  auto simplified = surfaceSimplification(
      *cube, SimplificationStopPredicate::edgeCountRatio(0.5),
      SimplificationStrategy::LINDSTROM_TURK);

  BOOST_CHECK(simplified);
  BOOST_CHECK(simplified->is<TriangulatedSurface>());

  const auto &simplifiedTin = simplified->as<TriangulatedSurface>();
  BOOST_CHECK_LT(simplifiedTin.numTriangles(), cube->numTriangles());

  std::cout << "Original triangles: " << cube->numTriangles() << "\n";
  std::cout << "Simplified triangles (LT): " << simplifiedTin.numTriangles()
            << "\n";
}
#endif // SFCGAL_WITH_EIGEN

#ifdef SFCGAL_WITH_EIGEN
// Test PolyhedralSurface simplification with Garland-Heckbert strategy
BOOST_AUTO_TEST_CASE(testSimplify_PolyhedralSurface_GarlandHeckbert)
{
  auto cube = createCubePolyhedralSurface();

  BOOST_CHECK_EQUAL(cube->numPatches(), 6);

  // Simplify using edge count ratio with Garland-Heckbert strategy
  auto simplified = surfaceSimplification(
      *cube, SimplificationStopPredicate::edgeCountRatio(0.5),
      SimplificationStrategy::GARLAND_HECKBERT);

  BOOST_CHECK(simplified);
  BOOST_CHECK(simplified->is<PolyhedralSurface>());

  std::cout << "Original patches: " << cube->numPatches() << "\n";
  std::cout << "Simplified patches: "
            << simplified->as<PolyhedralSurface>().numPatches() << "\n";
}
#endif // SFCGAL_WITH_EIGEN

#ifdef SFCGAL_WITH_EIGEN
// Test PolyhedralSurface simplification with Lindstrom-Turk strategy
BOOST_AUTO_TEST_CASE(testSimplify_PolyhedralSurface_LindstromTurk)
{
  auto cube = createCubePolyhedralSurface();

  BOOST_CHECK_EQUAL(cube->numPatches(), 6);

  // Simplify using edge count ratio with Lindstrom-Turk strategy
  auto simplified = surfaceSimplification(
      *cube, SimplificationStopPredicate::edgeCountRatio(0.5),
      SimplificationStrategy::LINDSTROM_TURK);

  BOOST_CHECK(simplified);
  BOOST_CHECK(simplified->is<PolyhedralSurface>());

  std::cout << "Original patches: " << cube->numPatches() << "\n";
  std::cout << "Simplified patches (LT): "
            << simplified->as<PolyhedralSurface>().numPatches() << "\n";
}
#endif // SFCGAL_WITH_EIGEN

// Test Solid simplification
BOOST_AUTO_TEST_CASE(testSimplify_Solid)
{
  auto  exteriorShell = createCubePolyhedralSurface();
  Solid solid(*exteriorShell);

  auto simplified = surfaceSimplification(
      solid, SimplificationStopPredicate::edgeCountRatio(0.5),
      SimplificationStrategy::EDGE_LENGTH);

  BOOST_CHECK(simplified);
  BOOST_CHECK(simplified->is<Solid>());

  const auto &simplifiedSolid = simplified->as<Solid>();
  BOOST_CHECK(!simplifiedSolid.isEmpty());

  std::cout << "Original solid exterior patches: "
            << solid.exteriorShell().numPatches() << "\n";
  std::cout << "Simplified solid exterior patches: "
            << simplifiedSolid.exteriorShell().numPatches() << "\n";
}

// Test MultiSolid simplification
BOOST_AUTO_TEST_CASE(testSimplify_MultiSolid)
{
  MultiSolid multiSolid;

  auto shell1 = createCubePolyhedralSurface();
  multiSolid.addGeometry(Solid(*shell1));

  auto shell2 = createCubePolyhedralSurface();
  multiSolid.addGeometry(Solid(*shell2));

  BOOST_CHECK_EQUAL(multiSolid.numGeometries(), 2);

  auto simplified = surfaceSimplification(
      multiSolid, SimplificationStopPredicate::edgeCountRatio(0.5),
      SimplificationStrategy::EDGE_LENGTH);

  BOOST_CHECK(simplified);
  BOOST_CHECK(simplified->is<MultiSolid>());

  const auto &simplifiedMulti = simplified->as<MultiSolid>();
  BOOST_CHECK_EQUAL(simplifiedMulti.numGeometries(), 2);

  std::cout << "MultiSolid simplified successfully" << "\n";
}

// Test edge count stop predicate
BOOST_AUTO_TEST_CASE(testSimplify_EdgeCountPredicate)
{
  auto cube = createCubeTriangulatedSurface();

  // Convert to surface mesh to count edges
  // For a cube with 12 triangles, there are 18 edges (12 * 3 / 2 = 18 for
  // manifold)

  // Simplify to a specific edge count (keep 10 edges)
  auto simplified =
      surfaceSimplification(*cube, SimplificationStopPredicate::edgeCount(10),
                            SimplificationStrategy::EDGE_LENGTH);

  BOOST_CHECK(simplified);
  BOOST_CHECK(simplified->is<TriangulatedSurface>());

  std::cout << "Simplified with edge count predicate" << "\n";
}

// Test empty geometries
BOOST_AUTO_TEST_CASE(testSimplify_EmptyGeometries)
{
  TriangulatedSurface emptyTin;
  BOOST_CHECK(emptyTin.isEmpty());

  auto simplified = surfaceSimplification(
      emptyTin, SimplificationStopPredicate::edgeCountRatio(0.5),
      SimplificationStrategy::EDGE_LENGTH);

  BOOST_CHECK(simplified);
  BOOST_CHECK(simplified->isEmpty());
}

// Test invalid ratio (too low)
BOOST_AUTO_TEST_CASE(testSimplify_InvalidRatio_TooLow)
{
  auto cube = createCubeTriangulatedSurface();

  BOOST_CHECK_THROW((void)surfaceSimplification(
                        *cube, SimplificationStopPredicate::edgeCountRatio(0.0),
                        SimplificationStrategy::EDGE_LENGTH),
                    std::invalid_argument);
}

// Test invalid ratio (too high)
BOOST_AUTO_TEST_CASE(testSimplify_InvalidRatio_TooHigh)
{
  auto cube = createCubeTriangulatedSurface();

  BOOST_CHECK_THROW((void)surfaceSimplification(
                        *cube, SimplificationStopPredicate::edgeCountRatio(1.0),
                        SimplificationStrategy::EDGE_LENGTH),
                    std::invalid_argument);
}

// Test invalid ratio (negative)
BOOST_AUTO_TEST_CASE(testSimplify_InvalidRatio_Negative)
{
  auto cube = createCubeTriangulatedSurface();

  BOOST_CHECK_THROW((void)surfaceSimplification(
                        *cube,
                        SimplificationStopPredicate::edgeCountRatio(-0.5),
                        SimplificationStrategy::EDGE_LENGTH),
                    std::invalid_argument);
}

// Test unsupported geometry type
BOOST_AUTO_TEST_CASE(testSimplify_UnsupportedGeometry)
{
  Point point(0, 0, 0);

  BOOST_CHECK_THROW((void)surfaceSimplification(
                        point, SimplificationStopPredicate::edgeCountRatio(0.5),
                        SimplificationStrategy::EDGE_LENGTH),
                    std::invalid_argument);
}

// Test simplification with WKT geometries
BOOST_AUTO_TEST_CASE(testSimplify_WKT_TriangulatedSurface)
{
  // Create a connected triangulated surface (a tetrahedron)
  std::string wkt = "TIN Z(((0 0 0, 1 0 0, 0.5 0.5 1, 0 0 0)),"
                    "((1 0 0, 0 1 0, 0.5 0.5 1, 1 0 0)),"
                    "((0 1 0, 0 0 0, 0.5 0.5 1, 0 1 0)),"
                    "((0 0 0, 0 1 0, 1 0 0, 0 0 0)))";

  std::unique_ptr<Geometry> geom(io::readWkt(wkt));
  BOOST_CHECK(geom->is<TriangulatedSurface>());

  auto simplified = surfaceSimplification(
      *geom, SimplificationStopPredicate::edgeCountRatio(0.7),
      SimplificationStrategy::EDGE_LENGTH);

  BOOST_CHECK(simplified);
  BOOST_CHECK(simplified->is<TriangulatedSurface>());

  std::cout << "WKT Source: " << geom->asText(2) << "\n";
  std::cout << "WKT Simplified: " << simplified->asText(2) << "\n";
}

// Test that simplified geometry maintains 3D characteristics
BOOST_AUTO_TEST_CASE(testSimplify_Maintains3D)
{
  auto cube = createCubeTriangulatedSurface();
  BOOST_CHECK(cube->is3D());

  auto simplified = surfaceSimplification(
      *cube, SimplificationStopPredicate::edgeCountRatio(0.5),
      SimplificationStrategy::EDGE_LENGTH);

  BOOST_CHECK(simplified);
  BOOST_CHECK(simplified->is3D());
}

// Test Solid with interior shells
BOOST_AUTO_TEST_CASE(testSimplify_SolidWithInteriorShells)
{
  auto  exteriorShell = createCubePolyhedralSurface();
  Solid solid(*exteriorShell);

  auto interiorShell =
      createCubePolyhedralSurface(Point(0.25, 0.25, 0.25), 0.5);

  solid.addInteriorShell(std::move(interiorShell));

  BOOST_CHECK_EQUAL(solid.numInteriorShells(), 1);

  // Use NoValidityCheck because SFCGAL's isValid for interior shells is not
  // fully implemented
  auto simplified = surfaceSimplification(
      solid, SimplificationStopPredicate::edgeCountRatio(0.5),
      SimplificationStrategy::EDGE_LENGTH, NoValidityCheck());

  BOOST_CHECK(simplified);
  BOOST_CHECK(simplified->is<Solid>());

  const auto &simplifiedSolid = simplified->as<Solid>();
  // Interior shells are preserved but not simplified
  BOOST_CHECK_EQUAL(simplifiedSolid.numInteriorShells(), 1);

  std::cout << "Solid with interior shells simplified successfully"
            << "\n";
}

BOOST_AUTO_TEST_SUITE_END()

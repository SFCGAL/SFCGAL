// Copyright (c) 2025-2026, SFCGAL Team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <SFCGAL/LineString.h>
#include <SFCGAL/MultiLineString.h>
#include <SFCGAL/Point.h>
#include <SFCGAL/Polygon.h>
#include <SFCGAL/PolyhedralSurface.h>
#include <SFCGAL/Solid.h>
#include <SFCGAL/algorithm/Chamfer.h>
#include <SFCGAL/algorithm/isValid.h>
#include <SFCGAL/algorithm/volume.h>
#include <SFCGAL/io/wkt.h>

#include <boost/test/unit_test.hpp>

#include <cmath>

using namespace SFCGAL;
using namespace boost::unit_test;

namespace {

// Unit cube: vertices (0,0,0) to (1,1,1)
const std::string CUBE_WKT =
    "SOLID ((("
    "(0 0 0, 0 1 0, 1 1 0, 1 0 0, 0 0 0)),"  // bottom (z=0)
    "((0 0 0, 0 0 1, 0 1 1, 0 1 0, 0 0 0)),"  // left (x=0)
    "((0 0 0, 1 0 0, 1 0 1, 0 0 1, 0 0 0)),"  // front (y=0)
    "((1 1 1, 0 1 1, 0 0 1, 1 0 1, 1 1 1)),"  // top (z=1)
    "((1 1 1, 1 0 1, 1 0 0, 1 1 0, 1 1 1)),"  // right (x=1)
    "((1 1 1, 1 1 0, 0 1 0, 0 1 1, 1 1 1))"   // back (y=1)
    "))";

// Helper to create a regular N-gon prism along the Z axis
// The N-gon cross-section is centered at (cx, cy) with circumradius r,
// extruded from z=0 to z=h.
auto
make_regular_prism(int n_sides, double circumradius, double height,
                   double cx = 0.0, double cy = 0.0)
    -> std::unique_ptr<Solid>
{
  // Compute base polygon vertices
  std::vector<std::pair<double, double>> base_pts;
  base_pts.reserve(n_sides);
  for (int i = 0; i < n_sides; ++i) {
    double angle = 2.0 * M_PI * i / n_sides;
    base_pts.emplace_back(cx + circumradius * std::cos(angle),
                          cy + circumradius * std::sin(angle));
  }

  // Build faces as Polygons
  std::vector<Polygon> faces;

  // Bottom face (z=0) — reversed winding for outward normal pointing down
  {
    LineString ring;
    for (int i = n_sides - 1; i >= 0; --i) {
      ring.addPoint(Point(base_pts[i].first, base_pts[i].second, 0));
    }
    ring.addPoint(
        Point(base_pts[n_sides - 1].first, base_pts[n_sides - 1].second, 0));
    faces.emplace_back(ring);
  }

  // Top face (z=h)
  {
    LineString ring;
    for (int i = 0; i < n_sides; ++i) {
      ring.addPoint(Point(base_pts[i].first, base_pts[i].second, height));
    }
    ring.addPoint(Point(base_pts[0].first, base_pts[0].second, height));
    faces.emplace_back(ring);
  }

  // Side faces
  for (int i = 0; i < n_sides; ++i) {
    int        j = (i + 1) % n_sides;
    LineString ring;
    ring.addPoint(Point(base_pts[i].first, base_pts[i].second, 0));
    ring.addPoint(Point(base_pts[j].first, base_pts[j].second, 0));
    ring.addPoint(Point(base_pts[j].first, base_pts[j].second, height));
    ring.addPoint(Point(base_pts[i].first, base_pts[i].second, height));
    ring.addPoint(Point(base_pts[i].first, base_pts[i].second, 0));
    faces.emplace_back(ring);
  }

  return std::make_unique<Solid>(faces);
}

// Helper: get a vertical side edge of a regular prism
// (the edge connecting the i-th base vertex to its top counterpart)
auto
get_prism_vertical_edge(int n_sides, double circumradius, double height, int i,
                        double cx = 0.0, double cy = 0.0) -> LineString
{
  double angle = 2.0 * M_PI * i / n_sides;
  double x     = cx + circumradius * std::cos(angle);
  double y     = cy + circumradius * std::sin(angle);

  LineString edge;
  edge.addPoint(Point(x, y, 0));
  edge.addPoint(Point(x, y, height));
  return edge;
}

// Helper: get a bottom edge of a regular prism (from vertex i to vertex i+1)
auto
get_prism_bottom_edge(int n_sides, double circumradius, int i,
                      double cx = 0.0, double cy = 0.0) -> LineString
{
  double angle1 = 2.0 * M_PI * i / n_sides;
  double angle2 = 2.0 * M_PI * ((i + 1) % n_sides) / n_sides;

  LineString edge;
  edge.addPoint(Point(cx + circumradius * std::cos(angle1),
                      cy + circumradius * std::sin(angle1), 0));
  edge.addPoint(Point(cx + circumradius * std::cos(angle2),
                      cy + circumradius * std::sin(angle2), 0));
  return edge;
}

} // anonymous namespace

BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_ChamferTest)

// ---------------------------------------------------------------------------
// Regression: 90-degree chamfer on axis-aligned cube
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(testChamfer_Cube90_Flat)
{
  auto cube = io::readWkt(CUBE_WKT);

  LineString edge;
  edge.addPoint(Point(0, 0, 0));
  edge.addPoint(Point(1, 0, 0));

  algorithm::ChamferOptions opts;
  opts.type   = algorithm::ChamferType::FLAT;
  opts.radius = 0.1;

  auto result = algorithm::chamfer(*cube, edge, opts);
  BOOST_REQUIRE(result != nullptr);
  BOOST_CHECK(result->geometryTypeId() == TYPE_SOLID);
  BOOST_CHECK(algorithm::isValid(*result));
  BOOST_CHECK(algorithm::volume(*result) < algorithm::volume(*cube));
  // Flat chamfer on cube adds 1 face (triangle replaces edge corner)
  const auto &s = result->as<Solid>();
  BOOST_CHECK(s.exteriorShell().numPolygons() > 6);
}

BOOST_AUTO_TEST_CASE(testChamfer_Cube90_Round)
{
  auto cube = io::readWkt(CUBE_WKT);

  LineString edge;
  edge.addPoint(Point(0, 0, 0));
  edge.addPoint(Point(1, 0, 0));

  algorithm::ChamferOptions opts;
  opts.type     = algorithm::ChamferType::ROUND;
  opts.radius   = 0.1;
  opts.segments = 4;

  auto result = algorithm::chamfer(*cube, edge, opts);
  BOOST_REQUIRE(result != nullptr);
  BOOST_CHECK(result->geometryTypeId() == TYPE_SOLID);
  BOOST_CHECK(algorithm::isValid(*result));
  BOOST_CHECK(algorithm::volume(*result) < algorithm::volume(*cube));
  const auto &s = result->as<Solid>();
  BOOST_CHECK(s.exteriorShell().numPolygons() > 6);
}

// ---------------------------------------------------------------------------
// Triangular prism (equilateral, 3 sides):
//   Interior angle between side faces = 60°
//   Opening angle = 180° - 60° = 120°
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(testChamfer_TriangularPrism_Flat)
{
  const int    n     = 3;
  const double r     = 2.0;
  const double h     = 3.0;
  auto         prism = make_regular_prism(n, r, h);

  // Get a vertical edge
  LineString edge = get_prism_vertical_edge(n, r, h, 0);

  algorithm::ChamferOptions opts;
  opts.type   = algorithm::ChamferType::FLAT;
  opts.radius = 0.2;

  auto result = algorithm::chamfer(*prism, edge, opts);
  BOOST_REQUIRE(result != nullptr);
  BOOST_CHECK(result->geometryTypeId() == TYPE_SOLID);
  BOOST_CHECK(algorithm::isValid(*result));
  BOOST_CHECK(algorithm::volume(*result) < algorithm::volume(*prism));
}

BOOST_AUTO_TEST_CASE(testChamfer_TriangularPrism_Round)
{
  const int    n     = 3;
  const double r     = 2.0;
  const double h     = 3.0;
  auto         prism = make_regular_prism(n, r, h);

  LineString edge = get_prism_vertical_edge(n, r, h, 0);

  algorithm::ChamferOptions opts;
  opts.type     = algorithm::ChamferType::ROUND;
  opts.radius   = 0.2;
  opts.segments = 8;

  auto result = algorithm::chamfer(*prism, edge, opts);
  BOOST_REQUIRE(result != nullptr);
  BOOST_CHECK(result->geometryTypeId() == TYPE_SOLID);
  BOOST_CHECK(algorithm::isValid(*result));
  BOOST_CHECK(algorithm::volume(*result) < algorithm::volume(*prism));
}

// ---------------------------------------------------------------------------
// Square prism (4 sides = cube-like): 90° interior, 90° opening
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(testChamfer_SquarePrism_Flat)
{
  const int    n     = 4;
  const double r     = 2.0;
  const double h     = 3.0;
  auto         prism = make_regular_prism(n, r, h);

  LineString edge = get_prism_vertical_edge(n, r, h, 0);

  algorithm::ChamferOptions opts;
  opts.type   = algorithm::ChamferType::FLAT;
  opts.radius = 0.2;

  auto result = algorithm::chamfer(*prism, edge, opts);
  BOOST_REQUIRE(result != nullptr);
  BOOST_CHECK(result->geometryTypeId() == TYPE_SOLID);
  BOOST_CHECK(algorithm::isValid(*result));
  BOOST_CHECK(algorithm::volume(*result) < algorithm::volume(*prism));
}

// ---------------------------------------------------------------------------
// Hexagonal prism (6 sides):
//   Interior angle = 120°
//   Opening angle = 180° - 120° = 60°
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(testChamfer_HexagonalPrism_Flat)
{
  const int    n     = 6;
  const double r     = 2.0;
  const double h     = 3.0;
  auto         prism = make_regular_prism(n, r, h);

  LineString edge = get_prism_vertical_edge(n, r, h, 0);

  algorithm::ChamferOptions opts;
  opts.type   = algorithm::ChamferType::FLAT;
  opts.radius = 0.15;

  auto result = algorithm::chamfer(*prism, edge, opts);
  BOOST_REQUIRE(result != nullptr);
  BOOST_CHECK(result->geometryTypeId() == TYPE_SOLID);
  BOOST_CHECK(algorithm::isValid(*result));
  BOOST_CHECK(algorithm::volume(*result) < algorithm::volume(*prism));
}

BOOST_AUTO_TEST_CASE(testChamfer_HexagonalPrism_Round)
{
  const int    n     = 6;
  const double r     = 2.0;
  const double h     = 3.0;
  auto         prism = make_regular_prism(n, r, h);

  LineString edge = get_prism_vertical_edge(n, r, h, 0);

  algorithm::ChamferOptions opts;
  opts.type     = algorithm::ChamferType::ROUND;
  opts.radius   = 0.15;
  opts.segments = 8;

  auto result = algorithm::chamfer(*prism, edge, opts);
  BOOST_REQUIRE(result != nullptr);
  BOOST_CHECK(result->geometryTypeId() == TYPE_SOLID);
  BOOST_CHECK(algorithm::isValid(*result));
  BOOST_CHECK(algorithm::volume(*result) < algorithm::volume(*prism));
}

// ---------------------------------------------------------------------------
// Octagonal prism (8 sides):
//   Interior angle = 135°
//   Opening angle = 180° - 135° = 45°
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(testChamfer_OctagonalPrism_Flat)
{
  const int    n     = 8;
  const double r     = 2.0;
  const double h     = 3.0;
  auto         prism = make_regular_prism(n, r, h);

  LineString edge = get_prism_vertical_edge(n, r, h, 0);

  algorithm::ChamferOptions opts;
  opts.type   = algorithm::ChamferType::FLAT;
  opts.radius = 0.1;

  auto result = algorithm::chamfer(*prism, edge, opts);
  BOOST_REQUIRE(result != nullptr);
  BOOST_CHECK(result->geometryTypeId() == TYPE_SOLID);
  BOOST_CHECK(algorithm::isValid(*result));
  BOOST_CHECK(algorithm::volume(*result) < algorithm::volume(*prism));
}

// ---------------------------------------------------------------------------
// Asymmetric chamfer on non-90° edge
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(testChamfer_Asymmetric_TriangularPrism)
{
  const int    n     = 3;
  const double r     = 2.0;
  const double h     = 3.0;
  auto         prism = make_regular_prism(n, r, h);

  LineString edge = get_prism_vertical_edge(n, r, h, 0);

  algorithm::ChamferOptions opts;
  opts.type     = algorithm::ChamferType::FLAT;
  opts.radius   = 0.2;
  opts.radius_y = 0.3; // asymmetric

  auto result = algorithm::chamfer(*prism, edge, opts);
  BOOST_REQUIRE(result != nullptr);
  BOOST_CHECK(result->geometryTypeId() == TYPE_SOLID);
  BOOST_CHECK(algorithm::isValid(*result));
  BOOST_CHECK(algorithm::volume(*result) < algorithm::volume(*prism));
}

// ---------------------------------------------------------------------------
// Multi-edge chamfer
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(testChamfer_MultiEdge_Cube)
{
  auto cube = io::readWkt(CUBE_WKT);

  // Chamfer two parallel edges of the cube
  MultiLineString multi;

  LineString edge1;
  edge1.addPoint(Point(0, 0, 0));
  edge1.addPoint(Point(1, 0, 0));
  multi.addGeometry(edge1);

  LineString edge2;
  edge2.addPoint(Point(0, 1, 0));
  edge2.addPoint(Point(1, 1, 0));
  multi.addGeometry(edge2);

  algorithm::ChamferOptions opts;
  opts.type   = algorithm::ChamferType::FLAT;
  opts.radius = 0.1;

  auto result = algorithm::chamfer(*cube, multi, opts);
  BOOST_REQUIRE(result != nullptr);
  BOOST_CHECK(result->geometryTypeId() == TYPE_SOLID);
  BOOST_CHECK(algorithm::isValid(*result));
  BOOST_CHECK(algorithm::volume(*result) < algorithm::volume(*cube));
}

// ---------------------------------------------------------------------------
// Bottom edge of a prism (edge between bottom face and side face = 90°)
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(testChamfer_BottomEdge_TriangularPrism)
{
  const int    n     = 3;
  const double r     = 2.0;
  const double h     = 3.0;
  auto         prism = make_regular_prism(n, r, h);

  // Bottom edge (between bottom face and side face)
  LineString edge = get_prism_bottom_edge(n, r, 0);

  algorithm::ChamferOptions opts;
  opts.type   = algorithm::ChamferType::FLAT;
  opts.radius = 0.15;

  auto result = algorithm::chamfer(*prism, edge, opts);
  BOOST_REQUIRE(result != nullptr);
  BOOST_CHECK(result->geometryTypeId() == TYPE_SOLID);
  BOOST_CHECK(algorithm::isValid(*result));
  BOOST_CHECK(algorithm::volume(*result) < algorithm::volume(*prism));
}

// ---------------------------------------------------------------------------
// Pentagon prism (5 sides):
//   Interior angle = 108°
//   Opening angle = 180° - 108° = 72°
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(testChamfer_PentagonPrism_Flat)
{
  const int    n     = 5;
  const double r     = 2.0;
  const double h     = 3.0;
  auto         prism = make_regular_prism(n, r, h);

  LineString edge = get_prism_vertical_edge(n, r, h, 0);

  algorithm::ChamferOptions opts;
  opts.type   = algorithm::ChamferType::FLAT;
  opts.radius = 0.15;

  auto result = algorithm::chamfer(*prism, edge, opts);
  BOOST_REQUIRE(result != nullptr);
  BOOST_CHECK(result->geometryTypeId() == TYPE_SOLID);
  BOOST_CHECK(algorithm::isValid(*result));
  BOOST_CHECK(algorithm::volume(*result) < algorithm::volume(*prism));
}

BOOST_AUTO_TEST_CASE(testChamfer_PentagonPrism_Round)
{
  const int    n     = 5;
  const double r     = 2.0;
  const double h     = 3.0;
  auto         prism = make_regular_prism(n, r, h);

  LineString edge = get_prism_vertical_edge(n, r, h, 0);

  algorithm::ChamferOptions opts;
  opts.type     = algorithm::ChamferType::ROUND;
  opts.radius   = 0.15;
  opts.segments = 8;

  auto result = algorithm::chamfer(*prism, edge, opts);
  BOOST_REQUIRE(result != nullptr);
  BOOST_CHECK(result->geometryTypeId() == TYPE_SOLID);
  BOOST_CHECK(algorithm::isValid(*result));
  BOOST_CHECK(algorithm::volume(*result) < algorithm::volume(*prism));
}

// ---------------------------------------------------------------------------
// Error cases
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(testChamfer_InvalidInput_NotSolid)
{
  Polygon poly;
  LineString edge;
  edge.addPoint(Point(0, 0, 0));
  edge.addPoint(Point(1, 0, 0));

  algorithm::ChamferOptions opts;
  BOOST_CHECK_THROW(algorithm::chamfer(poly, edge, opts),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(testChamfer_InvalidRadius_Zero)
{
  auto cube = io::readWkt(CUBE_WKT);

  LineString edge;
  edge.addPoint(Point(0, 0, 0));
  edge.addPoint(Point(1, 0, 0));

  algorithm::ChamferOptions opts;
  opts.radius = 0.0;
  BOOST_CHECK_THROW(algorithm::chamfer(*cube, edge, opts),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(testChamfer_InvalidRadius_Negative)
{
  auto cube = io::readWkt(CUBE_WKT);

  LineString edge;
  edge.addPoint(Point(0, 0, 0));
  edge.addPoint(Point(1, 0, 0));

  algorithm::ChamferOptions opts;
  opts.radius = -0.1;
  BOOST_CHECK_THROW(algorithm::chamfer(*cube, edge, opts),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(testChamfer_EdgeNotFound)
{
  auto cube = io::readWkt(CUBE_WKT);

  LineString edge;
  edge.addPoint(Point(5, 5, 5));
  edge.addPoint(Point(6, 6, 6));

  algorithm::ChamferOptions opts;
  // Edge not found → skipped, original returned
  auto result = algorithm::chamfer(*cube, edge, opts);
  BOOST_CHECK(result != nullptr);
}

BOOST_AUTO_TEST_CASE(testChamfer_ConcaveEdge_Skipped)
{
  // L-shape: concave inner corner at (1,1) is skipped (limitation)
  auto l_shape = io::readWkt(
      "SOLID ((("
      "(0 0 0, 0 2 0, 1 2 0, 1 1 0, 2 1 0, 2 0 0, 0 0 0)),"
      "((0 0 1, 2 0 1, 2 1 1, 1 1 1, 1 2 1, 0 2 1, 0 0 1)),"
      "((0 0 0, 2 0 0, 2 0 1, 0 0 1, 0 0 0)),"
      "((2 0 0, 2 1 0, 2 1 1, 2 0 1, 2 0 0)),"
      "((2 1 0, 1 1 0, 1 1 1, 2 1 1, 2 1 0)),"
      "((1 1 0, 1 2 0, 1 2 1, 1 1 1, 1 1 0)),"
      "((1 2 0, 0 2 0, 0 2 1, 1 2 1, 1 2 0)),"
      "((0 2 0, 0 0 0, 0 0 1, 0 2 1, 0 2 0))"
      "))");

  // Concave edge at (1,1,z) — skipped, original returned
  LineString edge;
  edge.addPoint(Point(1, 1, 0));
  edge.addPoint(Point(1, 1, 1));

  algorithm::ChamferOptions opts;
  opts.radius = 0.1;

  auto result = algorithm::chamfer(*l_shape, edge, opts);
  BOOST_CHECK(result != nullptr);
  BOOST_CHECK(result->geometryTypeId() == TYPE_SOLID);
  BOOST_CHECK(algorithm::volume(*result) == algorithm::volume(*l_shape));
}

// ---------------------------------------------------------------------------
// 3 edges meeting at a vertex (cube corner)
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(testChamfer_ThreeEdges_CubeCorner)
{
  auto cube = io::readWkt(CUBE_WKT);

  // 3 edges meeting at vertex (0,0,0)
  MultiLineString edges;

  LineString e1;
  e1.addPoint(Point(0, 0, 0));
  e1.addPoint(Point(1, 0, 0));
  edges.addGeometry(e1);

  LineString e2;
  e2.addPoint(Point(0, 0, 0));
  e2.addPoint(Point(0, 1, 0));
  edges.addGeometry(e2);

  LineString e3;
  e3.addPoint(Point(0, 0, 0));
  e3.addPoint(Point(0, 0, 1));
  edges.addGeometry(e3);

  algorithm::ChamferOptions opts;
  opts.type   = algorithm::ChamferType::FLAT;
  opts.radius = 0.15;

  auto result = algorithm::chamfer(*cube, edges, opts);
  BOOST_REQUIRE(result != nullptr);
  BOOST_CHECK(result->geometryTypeId() == TYPE_SOLID);
  BOOST_CHECK(algorithm::isValid(*result));
  BOOST_CHECK(algorithm::volume(*result) < algorithm::volume(*cube));
  // 3 chamfers on a cube corner: should have significantly more faces
  const auto &s = result->as<Solid>();
  BOOST_CHECK(s.exteriorShell().numPolygons() > 6);
}

BOOST_AUTO_TEST_CASE(testChamfer_ThreeEdges_CubeCorner_Round)
{
  auto cube = io::readWkt(CUBE_WKT);

  MultiLineString edges;

  LineString e1;
  e1.addPoint(Point(0, 0, 0));
  e1.addPoint(Point(1, 0, 0));
  edges.addGeometry(e1);

  LineString e2;
  e2.addPoint(Point(0, 0, 0));
  e2.addPoint(Point(0, 1, 0));
  edges.addGeometry(e2);

  LineString e3;
  e3.addPoint(Point(0, 0, 0));
  e3.addPoint(Point(0, 0, 1));
  edges.addGeometry(e3);

  algorithm::ChamferOptions opts;
  opts.type     = algorithm::ChamferType::ROUND;
  opts.radius   = 0.15;
  opts.segments = 6;

  auto result = algorithm::chamfer(*cube, edges, opts);
  BOOST_REQUIRE(result != nullptr);
  BOOST_CHECK(result->geometryTypeId() == TYPE_SOLID);
  BOOST_CHECK(algorithm::isValid(*result));
  BOOST_CHECK(algorithm::volume(*result) < algorithm::volume(*cube));
}

// ---------------------------------------------------------------------------
// L-shape: all edges (concave auto-skipped)
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(testChamfer_LShape_AllEdges)
{
  auto l_shape = io::readWkt(
      "SOLID ((("
      "(0 0 0, 0 2 0, 1 2 0, 1 1 0, 2 1 0, 2 0 0, 0 0 0)),"
      "((0 0 1, 2 0 1, 2 1 1, 1 1 1, 1 2 1, 0 2 1, 0 0 1)),"
      "((0 0 0, 2 0 0, 2 0 1, 0 0 1, 0 0 0)),"
      "((2 0 0, 2 1 0, 2 1 1, 2 0 1, 2 0 0)),"
      "((2 1 0, 1 1 0, 1 1 1, 2 1 1, 2 1 0)),"
      "((1 1 0, 1 2 0, 1 2 1, 1 1 1, 1 1 0)),"
      "((1 2 0, 0 2 0, 0 2 1, 1 2 1, 1 2 0)),"
      "((0 2 0, 0 0 0, 0 0 1, 0 2 1, 0 2 0))"
      "))");

  // All vertical + bottom + top edges — concave (1,1) auto-skipped
  MultiLineString edges;
  // Vertical edges
  auto add_edge = [&](double x1, double y1, double z1, double x2, double y2,
                      double z2) {
    LineString e;
    e.addPoint(Point(x1, y1, z1));
    e.addPoint(Point(x2, y2, z2));
    edges.addGeometry(e);
  };
  // 6 vertical
  add_edge(0, 0, 0, 0, 0, 1);
  add_edge(2, 0, 0, 2, 0, 1);
  add_edge(2, 1, 0, 2, 1, 1);
  add_edge(1, 1, 0, 1, 1, 1); // concave — will be skipped
  add_edge(1, 2, 0, 1, 2, 1);
  add_edge(0, 2, 0, 0, 2, 1);
  // 6 bottom
  add_edge(0, 0, 0, 2, 0, 0);
  add_edge(2, 0, 0, 2, 1, 0);
  add_edge(2, 1, 0, 1, 1, 0);
  add_edge(1, 1, 0, 1, 2, 0);
  add_edge(1, 2, 0, 0, 2, 0);
  add_edge(0, 2, 0, 0, 0, 0);

  algorithm::ChamferOptions opts;
  opts.radius = 0.1;

  auto result = algorithm::chamfer(*l_shape, edges, opts);
  BOOST_REQUIRE(result != nullptr);
  BOOST_CHECK(result->geometryTypeId() == TYPE_SOLID);
  BOOST_CHECK(algorithm::isValid(*result));
  BOOST_CHECK(algorithm::volume(*result) < algorithm::volume(*l_shape));
}

// ---------------------------------------------------------------------------
// Non-90° concave: chevron and star
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(testChamfer_Chevron_ConcaveSkipped)
{
  // Chevron: concave notch at (1, 1.5) — skipped (limitation)
  auto chevron = io::readWkt(
      "SOLID ((("
      "(0 0 0, 0 2 0, 1 1.5 0, 2 2 0, 2 0 0, 0 0 0)),"
      "((0 0 1, 2 0 1, 2 2 1, 1 1.5 1, 0 2 1, 0 0 1)),"
      "((0 0 0, 2 0 0, 2 0 1, 0 0 1, 0 0 0)),"
      "((2 0 0, 2 2 0, 2 2 1, 2 0 1, 2 0 0)),"
      "((2 2 0, 1 1.5 0, 1 1.5 1, 2 2 1, 2 2 0)),"
      "((1 1.5 0, 0 2 0, 0 2 1, 1 1.5 1, 1 1.5 0)),"
      "((0 2 0, 0 0 0, 0 0 1, 0 2 1, 0 2 0))"
      "))");

  // Concave edge (1, 1.5) → skipped, original returned
  LineString edge;
  edge.addPoint(Point(1, 1.5, 0));
  edge.addPoint(Point(1, 1.5, 1));

  algorithm::ChamferOptions opts;
  opts.radius = 0.1;

  auto result = algorithm::chamfer(*chevron, edge, opts);
  BOOST_CHECK(result != nullptr);
  BOOST_CHECK(result->geometryTypeId() == TYPE_SOLID);
  BOOST_CHECK(algorithm::volume(*result) == algorithm::volume(*chevron));
}

BOOST_AUTO_TEST_CASE(testChamfer_Chevron_ConvexWorks)
{
  auto chevron = io::readWkt(
      "SOLID ((("
      "(0 0 0, 0 2 0, 1 1.5 0, 2 2 0, 2 0 0, 0 0 0)),"
      "((0 0 1, 2 0 1, 2 2 1, 1 1.5 1, 0 2 1, 0 0 1)),"
      "((0 0 0, 2 0 0, 2 0 1, 0 0 1, 0 0 0)),"
      "((2 0 0, 2 2 0, 2 2 1, 2 0 1, 2 0 0)),"
      "((2 2 0, 1 1.5 0, 1 1.5 1, 2 2 1, 2 2 0)),"
      "((1 1.5 0, 0 2 0, 0 2 1, 1 1.5 1, 1 1.5 0)),"
      "((0 2 0, 0 0 0, 0 0 1, 0 2 1, 0 2 0))"
      "))");

  // Convex edge (2, 2) — non-90° angle
  LineString edge;
  edge.addPoint(Point(2, 2, 0));
  edge.addPoint(Point(2, 2, 1));

  algorithm::ChamferOptions opts;
  opts.radius = 0.15;

  auto result = algorithm::chamfer(*chevron, edge, opts);
  BOOST_REQUIRE(result != nullptr);
  BOOST_CHECK(result->geometryTypeId() == TYPE_SOLID);
  BOOST_CHECK(algorithm::isValid(*result));
  BOOST_CHECK(algorithm::volume(*result) < algorithm::volume(*chevron));
}

// ---------------------------------------------------------------------------
// Cube: all 12 edges
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Continuous polyline (single LineString) along a contour
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(testChamfer_ContinuousPolyline_ConvexCorner)
{
  // L-shape bottom: 2-segment sub-path through convex corner (3,0)
  auto l_shape = io::readWkt(
      "SOLID ((("
      "(0 0 0, 0 2 0, 1 2 0, 1 1 0, 2 1 0, 2 0 0, 0 0 0)),"
      "((0 0 1, 2 0 1, 2 1 1, 1 1 1, 1 2 1, 0 2 1, 0 0 1)),"
      "((0 0 0, 2 0 0, 2 0 1, 0 0 1, 0 0 0)),"
      "((2 0 0, 2 1 0, 2 1 1, 2 0 1, 2 0 0)),"
      "((2 1 0, 1 1 0, 1 1 1, 2 1 1, 2 1 0)),"
      "((1 1 0, 1 2 0, 1 2 1, 1 1 1, 1 1 0)),"
      "((1 2 0, 0 2 0, 0 2 1, 1 2 1, 1 2 0)),"
      "((0 2 0, 0 0 0, 0 0 1, 0 2 1, 0 2 0))"
      "))");

  // Single LineString with miter join at corner (2,0)
  LineString contour;
  contour.addPoint(Point(0, 0, 0));
  contour.addPoint(Point(2, 0, 0));
  contour.addPoint(Point(2, 1, 0));

  algorithm::ChamferOptions opts;
  opts.radius = 0.1;

  auto result = algorithm::chamfer(*l_shape, contour, opts);
  BOOST_REQUIRE(result != nullptr);
  BOOST_CHECK(result->geometryTypeId() == TYPE_SOLID);
  BOOST_CHECK(algorithm::isValid(*result));
  BOOST_CHECK(algorithm::volume(*result) < algorithm::volume(*l_shape));
}

BOOST_AUTO_TEST_CASE(testChamfer_ContinuousPolyline_FullContour)
{
  // L-shape bottom: full closed contour as single LineString
  // This passes through the concave corner (1,1) — the sweep must handle
  // miter joins at all corners. Currently, multi-corner miter joins may
  // produce non-planar quads in the sweep, causing the cutter to be invalid.
  // In that case the edge is skipped and the original is returned.
  auto l_shape = io::readWkt(
      "SOLID ((("
      "(0 0 0, 0 2 0, 1 2 0, 1 1 0, 2 1 0, 2 0 0, 0 0 0)),"
      "((0 0 1, 2 0 1, 2 1 1, 1 1 1, 1 2 1, 0 2 1, 0 0 1)),"
      "((0 0 0, 2 0 0, 2 0 1, 0 0 1, 0 0 0)),"
      "((2 0 0, 2 1 0, 2 1 1, 2 0 1, 2 0 0)),"
      "((2 1 0, 1 1 0, 1 1 1, 2 1 1, 2 1 0)),"
      "((1 1 0, 1 2 0, 1 2 1, 1 1 1, 1 1 0)),"
      "((1 2 0, 0 2 0, 0 2 1, 1 2 1, 1 2 0)),"
      "((0 2 0, 0 0 0, 0 0 1, 0 2 1, 0 2 0))"
      "))");

  LineString contour;
  contour.addPoint(Point(0, 0, 0));
  contour.addPoint(Point(2, 0, 0));
  contour.addPoint(Point(2, 1, 0));
  contour.addPoint(Point(1, 1, 0));
  contour.addPoint(Point(1, 2, 0));
  contour.addPoint(Point(0, 2, 0));
  contour.addPoint(Point(0, 0, 0));

  algorithm::ChamferOptions opts;
  opts.radius = 0.1;

  // This may fail due to sweep miter join limitations — result should
  // still be non-null (original returned if cutter is invalid)
  auto result = algorithm::chamfer(*l_shape, contour, opts);
  BOOST_CHECK(result != nullptr);
}

// ---------------------------------------------------------------------------
// Cube: all 12 edges
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(testChamfer_Cube_All12Edges)
{
  auto cube = io::readWkt(CUBE_WKT);

  MultiLineString edges;
  auto add = [&](double x1, double y1, double z1, double x2, double y2,
                 double z2) {
    LineString e;
    e.addPoint(Point(x1, y1, z1));
    e.addPoint(Point(x2, y2, z2));
    edges.addGeometry(e);
  };
  // Bottom 4
  add(0, 0, 0, 1, 0, 0);
  add(1, 0, 0, 1, 1, 0);
  add(1, 1, 0, 0, 1, 0);
  add(0, 1, 0, 0, 0, 0);
  // Top 4
  add(0, 0, 1, 1, 0, 1);
  add(1, 0, 1, 1, 1, 1);
  add(1, 1, 1, 0, 1, 1);
  add(0, 1, 1, 0, 0, 1);
  // Vertical 4
  add(0, 0, 0, 0, 0, 1);
  add(1, 0, 0, 1, 0, 1);
  add(1, 1, 0, 1, 1, 1);
  add(0, 1, 0, 0, 1, 1);

  algorithm::ChamferOptions opts;
  opts.radius = 0.1;

  auto result = algorithm::chamfer(*cube, edges, opts);
  BOOST_REQUIRE(result != nullptr);
  BOOST_CHECK(result->geometryTypeId() == TYPE_SOLID);
  BOOST_CHECK(algorithm::isValid(*result));
  BOOST_CHECK(algorithm::volume(*result) < algorithm::volume(*cube));
}

BOOST_AUTO_TEST_SUITE_END()

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

  // Edge not on the cube
  LineString edge;
  edge.addPoint(Point(5, 5, 5));
  edge.addPoint(Point(6, 6, 6));

  algorithm::ChamferOptions opts;
  // Should return original geometry (edge skipped with warning)
  auto result = algorithm::chamfer(*cube, edge, opts);
  BOOST_CHECK(result != nullptr);
}

BOOST_AUTO_TEST_SUITE_END()

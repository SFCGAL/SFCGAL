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

#include <CGAL/number_utils.h>

#include <boost/test/unit_test.hpp>

#include <cmath>
#include <vector>

using namespace SFCGAL;
using namespace boost::unit_test;

namespace {

// ─────────────────────────────────────────────────────────────────────────────
// Reference solids
// ─────────────────────────────────────────────────────────────────────────────

// Unit cube: vertices (0,0,0) to (1,1,1). Volume = 1.0.
// All edges are axis-aligned with 90° dihedral angle.
const std::string CUBE_WKT =
    "SOLID ((("
    "(0 0 0, 0 1 0, 1 1 0, 1 0 0, 0 0 0)),"  // bottom  (z=0, normal -Z)
    "((0 0 0, 0 0 1, 0 1 1, 0 1 0, 0 0 0))," // left    (x=0, normal -X)
    "((0 0 0, 1 0 0, 1 0 1, 0 0 1, 0 0 0))," // front   (y=0, normal -Y)
    "((1 1 1, 0 1 1, 0 0 1, 1 0 1, 1 1 1))," // top     (z=1, normal +Z)
    "((1 1 1, 1 0 1, 1 0 0, 1 1 0, 1 1 1))," // right   (x=1, normal +X)
    "((1 1 1, 1 1 0, 0 1 0, 0 1 1, 1 1 1))"  // back    (y=1, normal +Y)
    "))";

// 2×2×1 L-shape with concave inner corner at (1,1).
// Polygon cross-section: (0,0)→(2,0)→(2,1)→(1,1)→(1,2)→(0,2).
// Convex edges: (2,0,z), (0,0,z), (0,2,z), (1,2,z) (vertical);
//               bottom/top edges along convex outline.
// Concave edge: (1,1,z) (reflex corner → chamfer skipped automatically).
const std::string L_SHAPE_WKT =
    "SOLID ((("
    "(0 0 0, 0 2 0, 1 2 0, 1 1 0, 2 1 0, 2 0 0, 0 0 0))," // bottom
    "((0 0 1, 2 0 1, 2 1 1, 1 1 1, 1 2 1, 0 2 1, 0 0 1))," // top
    "((0 0 0, 2 0 0, 2 0 1, 0 0 1, 0 0 0)),"                // front (y=0)
    "((2 0 0, 2 1 0, 2 1 1, 2 0 1, 2 0 0)),"                // right (x=2)
    "((2 1 0, 1 1 0, 1 1 1, 2 1 1, 2 1 0)),"                // inner right (y=1)
    "((1 1 0, 1 2 0, 1 2 1, 1 1 1, 1 1 0)),"                // inner top (x=1)
    "((1 2 0, 0 2 0, 0 2 1, 1 2 1, 1 2 0)),"                // back (y=2)
    "((0 2 0, 0 0 0, 0 0 1, 0 2 1, 0 2 0))"                 // left (x=0)
    "))";

// Chevron prism: polygon (0,0)→(2,0)→(2,2)→(1,1.5)→(0,2), height=1.
// Convex vertical edges: (0,0,z), (2,0,z), (2,2,z), (0,2,z).
// Concave vertical edge: (1,1.5,z) — chamfer skipped automatically.
const std::string CHEVRON_WKT =
    "SOLID ((("
    "(0 0 0, 0 2 0, 1 1.5 0, 2 2 0, 2 0 0, 0 0 0))," // bottom
    "((0 0 1, 2 0 1, 2 2 1, 1 1.5 1, 0 2 1, 0 0 1))," // top
    "((0 0 0, 2 0 0, 2 0 1, 0 0 1, 0 0 0)),"            // front (y=0)
    "((2 0 0, 2 2 0, 2 2 1, 2 0 1, 2 0 0)),"            // right (x=2)
    "((2 2 0, 1 1.5 0, 1 1.5 1, 2 2 1, 2 2 0)),"        // notch right face
    "((1 1.5 0, 0 2 0, 0 2 1, 1 1.5 1, 1 1.5 0)),"      // notch left face
    "((0 2 0, 0 0 0, 0 0 1, 0 2 1, 0 2 0))"             // left (x=0)
    "))";

// ─────────────────────────────────────────────────────────────────────────────
// Helper: build a closed prism by extruding a 2D polygon to height h.
//
// @p pts  Polygon vertices in CCW order when viewed from above (+Z).
// @p h    Extrusion height.
//
// Face winding follows the SFCGAL convention: outward normals.
//   - Bottom (z=0): reversed order → normal pointing −Z.
//   - Top    (z=h): direct order  → normal pointing +Z.
//   - Sides       : (v_i, v_{i+1}, v_{i+1}+h, v_i+h) → outward normal.
// ─────────────────────────────────────────────────────────────────────────────
auto
make_prism(const std::vector<std::pair<double, double>> &pts, double h)
    -> std::unique_ptr<Solid>
{
  const int n = static_cast<int>(pts.size());

  std::vector<Polygon> faces;
  faces.reserve(n + 2);

  // Bottom face: reversed winding (normal pointing −Z)
  {
    LineString ring;
    for (int i = n - 1; i >= 0; --i) {
      ring.addPoint(Point(pts[i].first, pts[i].second, 0));
    }
    ring.addPoint(Point(pts[n - 1].first, pts[n - 1].second, 0));
    faces.emplace_back(ring);
  }

  // Top face: direct winding (normal pointing +Z)
  {
    LineString ring;
    for (int i = 0; i < n; ++i) {
      ring.addPoint(Point(pts[i].first, pts[i].second, h));
    }
    ring.addPoint(Point(pts[0].first, pts[0].second, h));
    faces.emplace_back(ring);
  }

  // Side faces: (v_i, v_{i+1}) bottom to top
  for (int i = 0; i < n; ++i) {
    const int j = (i + 1) % n;
    LineString ring;
    ring.addPoint(Point(pts[i].first, pts[i].second, 0));
    ring.addPoint(Point(pts[j].first, pts[j].second, 0));
    ring.addPoint(Point(pts[j].first, pts[j].second, h));
    ring.addPoint(Point(pts[i].first, pts[i].second, h));
    ring.addPoint(Point(pts[i].first, pts[i].second, 0));
    faces.emplace_back(ring);
  }

  return std::make_unique<Solid>(faces);
}

// Convenience: volume as double (avoids Kernel::FT in comparisons)
double
vol(const Geometry &g)
{
  return CGAL::to_double(algorithm::volume(g));
}

} // anonymous namespace

BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_ChamferTest)

// ═══════════════════════════════════════════════════════════════════════════
// SECTION A: Unit cube — 90° reference case
//
// All edges are axis-aligned with 90° dihedral angles. These tests establish
// the baseline behaviour of the chamfer algorithm.
// ═══════════════════════════════════════════════════════════════════════════

/**
 * Flat chamfer on one bottom edge of the unit cube.
 * Mirrors script example A01 (cube, radius=0.1).
 *
 * For a 90° corner with r=0.1 the removed volume is ≈ 0.5×0.1×0.1×1 = 0.005,
 * so the result volume should be just below 1.0 and above 0.98.
 */
BOOST_AUTO_TEST_CASE(testChamfer_Cube_Flat_BottomEdge)
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

  const double result_vol = vol(*result);
  BOOST_CHECK_LT(result_vol, vol(*cube)); // chamfer removed material
  BOOST_CHECK_GT(result_vol, 0.98);       // but not much (r=0.1 on unit cube)
}

/**
 * Round fillet on the same bottom edge.
 * Mirrors script example A02 (cube, type=1, radius=0.1, segments=8).
 */
BOOST_AUTO_TEST_CASE(testChamfer_Cube_Round_BottomEdge)
{
  auto cube = io::readWkt(CUBE_WKT);

  LineString edge;
  edge.addPoint(Point(0, 0, 0));
  edge.addPoint(Point(1, 0, 0));

  algorithm::ChamferOptions opts;
  opts.type     = algorithm::ChamferType::ROUND;
  opts.radius   = 0.1;
  opts.segments = 8;

  auto result = algorithm::chamfer(*cube, edge, opts);

  BOOST_REQUIRE(result != nullptr);
  BOOST_CHECK(result->geometryTypeId() == TYPE_SOLID);
  BOOST_CHECK(algorithm::isValid(*result));
  BOOST_CHECK_LT(vol(*result), vol(*cube));
}

/**
 * Flat chamfer on a vertical edge (Z-axis direction).
 * Mirrors script example A04 (cube, vertical edge at (1,0,0)→(1,0,1)).
 *
 * The sweep frame for a vertical axis requires the reference_normal = n1 fix:
 * T × Z_up is degenerate, the fallback T × Y_up gives the wrong frame
 * orientation for non-axis-aligned face normals.
 */
BOOST_AUTO_TEST_CASE(testChamfer_Cube_VerticalEdge)
{
  auto cube = io::readWkt(CUBE_WKT);

  // Vertical edge at x=1, y=0 (shared by front and right faces)
  LineString edge;
  edge.addPoint(Point(1, 0, 0));
  edge.addPoint(Point(1, 0, 1));

  algorithm::ChamferOptions opts;
  opts.type   = algorithm::ChamferType::FLAT;
  opts.radius = 0.1;

  auto result = algorithm::chamfer(*cube, edge, opts);

  BOOST_REQUIRE(result != nullptr);
  BOOST_CHECK(result->geometryTypeId() == TYPE_SOLID);
  BOOST_CHECK(algorithm::isValid(*result));
  BOOST_CHECK_LT(vol(*result), vol(*cube));
}

/**
 * Flat chamfer on all 4 bottom edges simultaneously (MultiLineString).
 * Mirrors script example A03 (cube, 4 bottom edges, radius=0.1).
 */
BOOST_AUTO_TEST_CASE(testChamfer_Cube_FourBottomEdges)
{
  auto cube = io::readWkt(CUBE_WKT);

  auto edges = io::readWkt(
      "MULTILINESTRING Z("
      "(0 0 0,1 0 0),(1 0 0,1 1 0),(1 1 0,0 1 0),(0 1 0,0 0 0))");

  algorithm::ChamferOptions opts;
  opts.type   = algorithm::ChamferType::FLAT;
  opts.radius = 0.1;

  auto result = algorithm::chamfer(*cube, *edges, opts);

  BOOST_REQUIRE(result != nullptr);
  BOOST_CHECK(result->geometryTypeId() == TYPE_SOLID);
  BOOST_CHECK(algorithm::isValid(*result));
  BOOST_CHECK_LT(vol(*result), vol(*cube));
}

/**
 * Three edges meeting at corner (0,0,0): flat and round chamfer.
 * Mirrors script examples in section F (3-edge corner chamfer on cube).
 *
 * All three adjacent cutters are unioned then subtracted from the cube,
 * testing the multi-cutter combination logic.
 */
BOOST_AUTO_TEST_CASE(testChamfer_Cube_ThreeEdgesAtCorner)
{
  auto edges = io::readWkt(
      "MULTILINESTRING Z("
      "(0 0 0,1 0 0),(0 0 0,0 1 0),(0 0 0,0 0 1))");

  // Flat
  {
    auto cube = io::readWkt(CUBE_WKT);
    algorithm::ChamferOptions opts;
    opts.type   = algorithm::ChamferType::FLAT;
    opts.radius = 0.1;

    auto result = algorithm::chamfer(*cube, *edges, opts);
    BOOST_REQUIRE(result != nullptr);
    BOOST_CHECK(result->geometryTypeId() == TYPE_SOLID);
    BOOST_CHECK(algorithm::isValid(*result));
    BOOST_CHECK_LT(vol(*result), vol(*cube));
  }

  // Round
  {
    auto cube = io::readWkt(CUBE_WKT);
    algorithm::ChamferOptions opts;
    opts.type     = algorithm::ChamferType::ROUND;
    opts.radius   = 0.1;
    opts.segments = 6;

    auto result = algorithm::chamfer(*cube, *edges, opts);
    BOOST_REQUIRE(result != nullptr);
    BOOST_CHECK(result->geometryTypeId() == TYPE_SOLID);
    BOOST_CHECK(algorithm::isValid(*result));
    BOOST_CHECK_LT(vol(*result), vol(*cube));
  }
}

// ═══════════════════════════════════════════════════════════════════════════
// SECTION B: Non-axis-aligned geometry
//
// The chamfer sweep frame is computed as T × Z_up. When solid faces are
// rotated relative to the coordinate axes the resulting frame normal can be
// anti-parallel to n1, placing the profile outside the solid. The fix
// reference_normal = n1 corrects the frame orientation for single-segment edges.
// ═══════════════════════════════════════════════════════════════════════════

/**
 * Flat chamfer on a vertical edge of a diamond prism (45°-rotated square).
 * Polygon cross-section: (2,0)→(0,2)→(−2,0)→(0,−2), height=3.
 * Mirrors script example D section (diamond_flat).
 *
 * The dihedral angle is 90° (same as the cube) but the two adjacent faces are
 * rotated 45° relative to XY. Without reference_normal = n1, the chamfer
 * profile is placed on the wrong side and the boolean difference is near-zero.
 */
BOOST_AUTO_TEST_CASE(testChamfer_DiamondPrism_Flat)
{
  auto prism = make_prism({{2, 0}, {0, 2}, {-2, 0}, {0, -2}}, 3.0);

  const double vol_original = vol(*prism);

  // Vertical edge at the rightmost vertex (2,0)
  LineString edge;
  edge.addPoint(Point(2, 0, 0));
  edge.addPoint(Point(2, 0, 3));

  algorithm::ChamferOptions opts;
  opts.type   = algorithm::ChamferType::FLAT;
  opts.radius = 0.2;

  auto result = algorithm::chamfer(*prism, edge, opts);

  BOOST_REQUIRE(result != nullptr);
  BOOST_CHECK(result->geometryTypeId() == TYPE_SOLID);
  BOOST_CHECK(algorithm::isValid(*result));
  BOOST_CHECK_LT(vol(*result), vol_original);
}

// ═══════════════════════════════════════════════════════════════════════════
// SECTION C: Non-90° dihedral angles
// ═══════════════════════════════════════════════════════════════════════════

/**
 * Flat chamfer on a vertical edge of an equilateral triangular prism.
 * 3-sided polygon, circumradius=2, height=3.
 * Interior angle 60° → dihedral opening angle ≈ 120°.
 * Mirrors script example B section (3-gon flat chamfer).
 */
BOOST_AUTO_TEST_CASE(testChamfer_TriangularPrism_Flat)
{
  const double r = 2.0;
  auto         prism =
      make_prism({{r, 0},
                  {r * std::cos(2 * M_PI / 3), r * std::sin(2 * M_PI / 3)},
                  {r * std::cos(4 * M_PI / 3), r * std::sin(4 * M_PI / 3)}},
                 3.0);

  const double vol_original = vol(*prism);

  LineString edge;
  edge.addPoint(Point(r, 0, 0));
  edge.addPoint(Point(r, 0, 3));

  algorithm::ChamferOptions opts;
  opts.type   = algorithm::ChamferType::FLAT;
  opts.radius = 0.2;

  auto result = algorithm::chamfer(*prism, edge, opts);

  BOOST_REQUIRE(result != nullptr);
  BOOST_CHECK(result->geometryTypeId() == TYPE_SOLID);
  BOOST_CHECK(algorithm::isValid(*result));
  BOOST_CHECK_LT(vol(*result), vol_original);
}

/**
 * Flat chamfer on a vertical edge of a regular hexagonal prism.
 * 6 sides, circumradius=2, height=3.
 * Interior angle 120° → dihedral opening angle ≈ 60°.
 * Mirrors script example B section (6-gon flat chamfer).
 */
BOOST_AUTO_TEST_CASE(testChamfer_HexagonalPrism_Flat)
{
  const double r = 2.0;
  const int    n = 6;

  std::vector<std::pair<double, double>> pts(n);
  for (int i = 0; i < n; ++i) {
    pts[i] = {r * std::cos(2 * M_PI * i / n), r * std::sin(2 * M_PI * i / n)};
  }
  auto prism = make_prism(pts, 3.0);

  const double vol_original = vol(*prism);

  LineString edge;
  edge.addPoint(Point(r, 0, 0));
  edge.addPoint(Point(r, 0, 3));

  algorithm::ChamferOptions opts;
  opts.type   = algorithm::ChamferType::FLAT;
  opts.radius = 0.15;

  auto result = algorithm::chamfer(*prism, edge, opts);

  BOOST_REQUIRE(result != nullptr);
  BOOST_CHECK(result->geometryTypeId() == TYPE_SOLID);
  BOOST_CHECK(algorithm::isValid(*result));
  BOOST_CHECK_LT(vol(*result), vol_original);
}

// ═══════════════════════════════════════════════════════════════════════════
// SECTION D: L-shape — concave auto-skip and multi-edge scenarios
//
// The L-shape has a reflex (concave) corner at (1,1). The chamfer algorithm
// detects and skips concave edges automatically, always returning a valid solid.
// ═══════════════════════════════════════════════════════════════════════════

/**
 * Chamfer on the concave (reflex) edge at (1,1) of the L-shape.
 * The edge is detected as concave and skipped → result is identical to input.
 */
BOOST_AUTO_TEST_CASE(testChamfer_LShape_ConcaveEdge_Skipped)
{
  auto l_shape = io::readWkt(L_SHAPE_WKT);

  LineString edge; // concave vertical edge at the inner corner
  edge.addPoint(Point(1, 1, 0));
  edge.addPoint(Point(1, 1, 1));

  algorithm::ChamferOptions opts;
  opts.radius = 0.1;

  auto result = algorithm::chamfer(*l_shape, edge, opts);

  BOOST_REQUIRE(result != nullptr);
  BOOST_CHECK(result->geometryTypeId() == TYPE_SOLID);
  // Skipped edge → volume unchanged (clone returned)
  BOOST_CHECK_CLOSE(vol(*result), vol(*l_shape), 1e-10);
}

/**
 * Chamfer all 18 edges of the L-shape (6 vertical + 6 bottom + 6 top).
 * The concave edges at (1,1) are auto-skipped; the remaining 14 are chamfered.
 * Mirrors script examples G section (L-shape all edges).
 */
BOOST_AUTO_TEST_CASE(testChamfer_LShape_AllEdges)
{
  auto l_shape = io::readWkt(L_SHAPE_WKT);

  // All 18 edges — the 4 at (1,1) will be silently skipped as concave
  auto edges = io::readWkt(
      "MULTILINESTRING Z("
      "(0 0 0,0 0 1),(2 0 0,2 0 1),(2 1 0,2 1 1),"
      "(1 1 0,1 1 1),"            // concave — skipped
      "(1 2 0,1 2 1),(0 2 0,0 2 1),"
      "(0 0 0,2 0 0),(2 0 0,2 1 0),(2 1 0,1 1 0),"
      "(1 1 0,1 2 0),"            // concave — skipped
      "(1 2 0,0 2 0),(0 2 0,0 0 0),"
      "(0 0 1,2 0 1),(2 0 1,2 1 1),(2 1 1,1 1 1),"
      "(1 1 1,1 2 1),"            // concave — skipped
      "(1 2 1,0 2 1),(0 2 1,0 0 1))");

  algorithm::ChamferOptions opts;
  opts.type   = algorithm::ChamferType::FLAT;
  opts.radius = 0.08;

  auto result = algorithm::chamfer(*l_shape, *edges, opts);

  BOOST_REQUIRE(result != nullptr);
  BOOST_CHECK(result->geometryTypeId() == TYPE_SOLID);
  BOOST_CHECK(algorithm::isValid(*result));
  BOOST_CHECK_LT(vol(*result), vol(*l_shape));
}

/**
 * Chamfer via a multi-segment LineString (convex sub-contour with miter join).
 * The LineString Z(0 0 0, 2 0 0, 2 1 0) covers two bottom edges of the L-shape,
 * joined at the convex corner (2,0,0) with a bisector-plane miter join.
 * Mirrors script example G section (L-shape contour convex, miter join).
 */
BOOST_AUTO_TEST_CASE(testChamfer_LShape_ContourMiterJoin)
{
  auto l_shape = io::readWkt(L_SHAPE_WKT);

  // 2-segment bottom sub-contour with miter join at corner (2,0,0)
  LineString contour;
  contour.addPoint(Point(0, 0, 0));
  contour.addPoint(Point(2, 0, 0));
  contour.addPoint(Point(2, 1, 0));

  algorithm::ChamferOptions opts;
  opts.type   = algorithm::ChamferType::FLAT;
  opts.radius = 0.1;

  auto result = algorithm::chamfer(*l_shape, contour, opts);

  BOOST_REQUIRE(result != nullptr);
  BOOST_CHECK(result->geometryTypeId() == TYPE_SOLID);
  BOOST_CHECK(algorithm::isValid(*result));
  BOOST_CHECK_LT(vol(*result), vol(*l_shape));
}

// ═══════════════════════════════════════════════════════════════════════════
// SECTION E: Non-rectangular polygon (chevron)
// ═══════════════════════════════════════════════════════════════════════════

/**
 * Flat chamfer on a convex edge of the chevron prism.
 * The faces adjacent to vertex (2,2,z) form a non-90° dihedral angle.
 * Mirrors script example H section (chevron convex corner).
 */
BOOST_AUTO_TEST_CASE(testChamfer_Chevron_ConvexEdge)
{
  auto chevron = io::readWkt(CHEVRON_WKT);

  // Convex corner (2,2,z)
  LineString edge;
  edge.addPoint(Point(2, 2, 0));
  edge.addPoint(Point(2, 2, 1));

  algorithm::ChamferOptions opts;
  opts.type   = algorithm::ChamferType::FLAT;
  opts.radius = 0.15;

  auto result = algorithm::chamfer(*chevron, edge, opts);

  BOOST_REQUIRE(result != nullptr);
  BOOST_CHECK(result->geometryTypeId() == TYPE_SOLID);
  BOOST_CHECK(algorithm::isValid(*result));
  BOOST_CHECK_LT(vol(*result), vol(*chevron));
}

// ═══════════════════════════════════════════════════════════════════════════
// SECTION F: Error cases
// ═══════════════════════════════════════════════════════════════════════════

/**
 * chamfer() requires a Solid or PolyhedralSurface. Passing a Polygon throws.
 */
BOOST_AUTO_TEST_CASE(testChamfer_Error_NotSolid)
{
  Polygon    poly;
  LineString edge;
  edge.addPoint(Point(0, 0, 0));
  edge.addPoint(Point(1, 0, 0));

  algorithm::ChamferOptions opts;
  BOOST_CHECK_THROW(algorithm::chamfer(poly, edge, opts), std::invalid_argument);
}

/**
 * Radius ≤ 0 must throw std::invalid_argument.
 */
BOOST_AUTO_TEST_CASE(testChamfer_Error_InvalidRadius)
{
  auto cube = io::readWkt(CUBE_WKT);

  LineString edge;
  edge.addPoint(Point(0, 0, 0));
  edge.addPoint(Point(1, 0, 0));

  algorithm::ChamferOptions opts;

  opts.radius = 0.0;
  BOOST_CHECK_THROW(algorithm::chamfer(*cube, edge, opts), std::invalid_argument);

  opts.radius = -0.1;
  BOOST_CHECK_THROW(algorithm::chamfer(*cube, edge, opts), std::invalid_argument);
}

/**
 * Edge not found in the solid mesh → silently skipped, original returned.
 * Must not throw.
 */
BOOST_AUTO_TEST_CASE(testChamfer_Error_EdgeNotFound)
{
  auto cube = io::readWkt(CUBE_WKT);

  LineString edge;
  edge.addPoint(Point(5, 5, 5)); // far outside the unit cube
  edge.addPoint(Point(6, 6, 6));

  algorithm::ChamferOptions opts;
  opts.radius = 0.1;

  auto result = algorithm::chamfer(*cube, edge, opts);
  BOOST_CHECK(result != nullptr);
}

BOOST_AUTO_TEST_SUITE_END()

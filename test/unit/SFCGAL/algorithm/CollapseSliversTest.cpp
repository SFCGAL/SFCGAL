// Copyright (c) 2026-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <boost/test/unit_test.hpp>

#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiSolid.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/algorithm/Chamfer.h"
#include "SFCGAL/algorithm/difference.h"
#include "SFCGAL/algorithm/eliminateSlivers.h"
#include "SFCGAL/algorithm/isValid.h"
#include "SFCGAL/algorithm/volume.h"
#include "SFCGAL/io/wkt.h"
#include "SFCGAL/primitive3d/Cylinder.h"

#include <CGAL/number_utils.h>

#include <cmath>
#include <limits>

using namespace SFCGAL;
using namespace boost::unit_test;

namespace {

/// Shortest edge of every ring of every patch of every shell.
auto
shortestEdge(const Solid &solid) -> double
{
  double shortest = std::numeric_limits<double>::max();

  for (std::size_t s = 0; s < solid.numShells(); ++s) {
    const PolyhedralSurface &shell = solid.shellN(s);
    for (std::size_t p = 0; p < shell.numPatches(); ++p) {
      const Polygon &patch = shell.patchN(p);
      for (std::size_t r = 0; r < patch.numRings(); ++r) {
        const LineString &ring = patch.ringN(r);
        for (std::size_t i = 0; i + 1 < ring.numPoints(); ++i) {
          const Kernel::Vector_3 edge =
              ring.pointN(i + 1).toVector_3() - ring.pointN(i).toVector_3();
          shortest = std::min(
              shortest, std::sqrt(CGAL::to_double(edge.squared_length())));
        }
      }
    }
  }
  return shortest;
}

} // namespace

BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_CollapseSliversTest)

/**
 * The chamfer cutter is built from floating point trigonometry, so it only
 * touches the box it cuts approximately: the exact corefinement then emits
 * vertices a rounding error apart, and needle faces between them. This box has
 * integer coordinates, so nothing but that construction can introduce them.
 *
 * Without collapseSlivers() the result carries edges of about 2e-16, which is
 * valid on its own but degenerates as soon as the shell is triangulated again
 * or written with a coarser precision.
 */
BOOST_AUTO_TEST_CASE(testBooleanLeavesNoSliverEdge)
{
  std::unique_ptr<Geometry> const box(
      io::readWkt("SOLID Z ((((0 0 0,0 4 0,4 4 0,4 0 0,0 0 0)),"
                  "((0 0 3,4 0 3,4 4 3,0 4 3,0 0 3)),"
                  "((0 0 0,4 0 0,4 0 3,0 0 3,0 0 0)),"
                  "((4 0 0,4 4 0,4 4 3,4 0 3,4 0 0)),"
                  "((4 4 0,0 4 0,0 4 3,4 4 3,4 4 0)),"
                  "((0 4 0,0 0 0,0 0 3,0 4 3,0 4 0))))"));
  std::unique_ptr<Geometry> const edge(
      io::readWkt("LINESTRING Z (4 0 0,4 0 3)"));

  algorithm::ChamferOptions options;
  options.type   = algorithm::ChamferType::FLAT;
  options.radius = 0.5;

  std::unique_ptr<Geometry> const result =
      algorithm::chamfer(*box, *edge, options);
  BOOST_REQUIRE(result != nullptr);
  BOOST_REQUIRE(result->geometryTypeId() == TYPE_SOLID);

  // The chamfer face is 0.35 wide, everything else is of the order of the box
  BOOST_CHECK_GT(shortestEdge(result->as<Solid>()), 1e-6);

  BOOST_CHECK(algorithm::isValid(*result));
  BOOST_CHECK_LT(algorithm::volume(result->as<Solid>()), Kernel::FT(48));
}

/**
 * Rounding the coordinates before the operation breaks the exact incidences
 * the two cylinders had, and the corefinement then lands its intersection
 * points a rounding away from the vertices they used to meet. Those slivers
 * are as large as the rounding, far above what the boolean cleans up on its
 * own, and only the caller knows how coarse the rounding was.
 */
BOOST_AUTO_TEST_CASE(testEliminateSliversTakesAnExplicitTolerance)
{
  const long                scale = 1000000; // six decimals
  std::unique_ptr<Geometry> outer(
      Cylinder(100.0, 30.0).generatePolyhedralSurface().clone());
  std::unique_ptr<Geometry> inner(
      Cylinder(70.0, 30.0).generatePolyhedralSurface().clone());
  outer->round(scale);
  inner->round(scale);

  const Solid outerSolid(outer->as<PolyhedralSurface>());
  const Solid innerSolid(inner->as<PolyhedralSurface>());

  std::unique_ptr<Geometry> const ring =
      algorithm::difference3D(outerSolid, innerSolid);
  BOOST_REQUIRE(ring->geometryTypeId() == TYPE_SOLID);

  // The rounding is 1e-6, the slivers it leaves are of that order
  const double before = shortestEdge(ring->as<Solid>());
  BOOST_CHECK_LT(before, 1e-6);

  // The default bound targets the rounding of the construction, not this one
  std::unique_ptr<Geometry> const untouched =
      algorithm::eliminateSlivers(*ring);
  BOOST_CHECK_LT(shortestEdge(untouched->as<Solid>()), 1e-6);

  std::unique_ptr<Geometry> const repaired =
      algorithm::eliminateSlivers(*ring, 1e-6);
  BOOST_CHECK_GT(shortestEdge(repaired->as<Solid>()), 1e-6);

  BOOST_CHECK(algorithm::isValid(*repaired));
  BOOST_CHECK_EQUAL(algorithm::volume(repaired->as<Solid>()),
                    algorithm::volume(ring->as<Solid>()));
}

/**
 * Writing a geometry with a precision coarser than its own slivers folds their
 * rings onto themselves: both ends of the needle land on the same coordinates.
 * Such a patch has no plane left, so it cannot even be triangulated, and
 * eliminateSlivers() has to drop it before anything else.
 */
BOOST_AUTO_TEST_CASE(testEliminateSliversUnfoldsARoundedResult)
{
  const long                scale = 1000000; // six decimals
  std::unique_ptr<Geometry> outer(
      Cylinder(100.0, 30.0).generatePolyhedralSurface().clone());
  std::unique_ptr<Geometry> inner(
      Cylinder(70.0, 30.0).generatePolyhedralSurface().clone());
  outer->round(scale);
  inner->round(scale);

  std::unique_ptr<Geometry> ring =
      algorithm::difference3D(Solid(outer->as<PolyhedralSurface>()),
                              Solid(inner->as<PolyhedralSurface>()));
  BOOST_REQUIRE(ring->geometryTypeId() == TYPE_SOLID);

  // Same rounding on the result: the slivers it holds are smaller than that
  ring->round(scale);
  BOOST_CHECK(!algorithm::isValid(*ring));

  std::unique_ptr<Geometry> const repaired = algorithm::eliminateSlivers(*ring);
  BOOST_CHECK(algorithm::isValid(*repaired));
  BOOST_CHECK_LT(repaired->as<Solid>().exteriorShell().numPatches(),
                 ring->as<Solid>().exteriorShell().numPatches());
}

/// A cap is a flat triangle, not a needle: only the negligible ones may be
/// flipped, so a wide one has to survive untouched.
BOOST_AUTO_TEST_CASE(testEliminateSliversKeepsAWideCap)
{
  // Tetrahedron whose base is a cap: the middle vertex sits 1 mm off the line
  // joining the two others, an angle of almost 180 degrees
  std::unique_ptr<Geometry> const tetrahedron(
      io::readWkt("SOLID Z ((((0 0 0,5 0.001 0,10 0 0,0 0 0)),"
                  "((0 0 0,10 0 0,5 0.0005 1,0 0 0)),"
                  "((10 0 0,5 0.001 0,5 0.0005 1,10 0 0)),"
                  "((5 0.001 0,0 0 0,5 0.0005 1,5 0.001 0))))"));
  BOOST_REQUIRE(algorithm::isValid(*tetrahedron));

  std::unique_ptr<Geometry> const repaired =
      algorithm::eliminateSlivers(*tetrahedron);

  BOOST_CHECK_EQUAL(repaired->as<Solid>().exteriorShell().numPatches(),
                    tetrahedron->as<Solid>().exteriorShell().numPatches());
  BOOST_CHECK_EQUAL(algorithm::volume(repaired->as<Solid>()),
                    algorithm::volume(tetrahedron->as<Solid>()));
  BOOST_CHECK(algorithm::isValid(*repaired));
}

/// A collection keeps its type, a MultiSolid does not become a plain one.
BOOST_AUTO_TEST_CASE(testEliminateSliversKeepsTheCollectionType)
{
  std::unique_ptr<Geometry> const multiSolid(io::readWkt(
      "MULTISOLID Z (((((0 0 0,0 1 0,1 1 0,1 0 0,0 0 0)),"
      "((0 0 1,1 0 1,1 1 1,0 1 1,0 0 1)),((0 0 0,1 0 0,1 0 1,0 0 1,0 0 0)),"
      "((1 0 0,1 1 0,1 1 1,1 0 1,1 0 0)),((1 1 0,0 1 0,0 1 1,1 1 1,1 1 0)),"
      "((0 1 0,0 0 0,0 0 1,0 1 1,0 1 0)))))"));
  BOOST_REQUIRE(multiSolid->geometryTypeId() == TYPE_MULTISOLID);

  std::unique_ptr<Geometry> const repaired =
      algorithm::eliminateSlivers(*multiSolid);

  BOOST_CHECK(repaired->geometryTypeId() == TYPE_MULTISOLID);
  BOOST_CHECK_EQUAL(repaired->numGeometries(), 1U);
  BOOST_CHECK(repaired->geometryN(0).geometryTypeId() == TYPE_SOLID);
}

/// Nothing to repair on a geometry that holds no face at all.
BOOST_AUTO_TEST_CASE(testEliminateSliversLeavesOtherTypesAlone)
{
  const std::string wkts[] = {"POINT (1 2 3)", "LINESTRING (0 0,1 1)",
                              "POLYGON ((0 0,1 0,1 1,0 1,0 0))",
                              "GEOMETRYCOLLECTION EMPTY"};

  for (const std::string &wkt : wkts) {
    std::unique_ptr<Geometry> const geometry(io::readWkt(wkt));
    std::unique_ptr<Geometry> const repaired =
        algorithm::eliminateSlivers(*geometry);
    BOOST_CHECK_MESSAGE(repaired->almostEqual(*geometry, -1),
                        wkt << " came back as " << repaired->asText(1));
  }
}

BOOST_AUTO_TEST_SUITE_END()

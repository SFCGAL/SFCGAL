// Copyright (c) 2026-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <boost/test/unit_test.hpp>

#include "SFCGAL/LineString.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/algorithm/Chamfer.h"
#include "SFCGAL/algorithm/isValid.h"
#include "SFCGAL/algorithm/volume.h"
#include "SFCGAL/io/wkt.h"

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

BOOST_AUTO_TEST_SUITE_END()

// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <boost/test/unit_test.hpp>

#include <array>

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/Kernel.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiLineString.h"
#include "SFCGAL/MultiPoint.h"
#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/MultiSolid.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/TriangulatedSurface.h"
#include "SFCGAL/algorithm/plane.h"
#include "SFCGAL/io/wkt.h"

using namespace SFCGAL;
using namespace boost::unit_test;

BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_PlaneTest)

BOOST_AUTO_TEST_CASE(testPlane1)
{
  std::unique_ptr<Geometry> gA(io::readWkt("POLYGON ((0 0,1 0,1 1,0 1,0 0))"));

  CGAL::Plane_3<Kernel> const plane =
      algorithm::plane3D<Kernel>(gA->as<Polygon>());
  BOOST_CHECK_EQUAL(plane.a(), 0.0);
  BOOST_CHECK_EQUAL(plane.b(), 0.0);
  BOOST_CHECK_EQUAL(plane.c(), 2.0);
}

BOOST_AUTO_TEST_CASE(testPlane)
{
  struct TestCase {
    const std::string _wkt;
    const bool        _isPlane;
  };
  const std::array<TestCase, 11> test = {{
      // only two points
      {._wkt = "LINESTRING (1 2 3,4 5 6)", ._isPlane = true},
      // all points in the same place
      {._wkt = "LINESTRING (1 2 3,1 2 3,1 2 3,1 2 3)", ._isPlane = true},
      // all points aligned
      {._wkt = "LINESTRING (1 2 3,2 4 6,3 6 9,4 8 12)", ._isPlane = true},
      // triangle must be plane
      {._wkt = "LINESTRING (1 2 3,6 5 4,7 8 9)", ._isPlane = true},
      // all point in the plane z=0
      {._wkt = "LINESTRING (0 0 0,1 0 0,1 1 0,0 1 0,0 0 0)", ._isPlane = true},
      // all points in the plane x=2
      {._wkt = "LINESTRING (2 1 0,2 0 0,2 1 0,2 1 0,2 0 3)", ._isPlane = true},
      // one point out of plane
      {._wkt = "LINESTRING (2 1 0,2 0 0,2 1 1,2 1 0,1 0 3)", ._isPlane = false},
      // fix #247
      {._wkt     = "LINESTRING (0 0 0, 1e-5 0 0, 1e-5 1e-5 0, 0 1e-5 1e-5)",
       ._isPlane = false},
      // self crossing outline: Newell's sum vanishes, the three point estimate
      // takes over
      {._wkt = "LINESTRING (0 0 0,1 1 0,1 0 0,0 1 0,0 0 0)", ._isPlane = true},
      // tilted plane x = z, with three nearly collinear vertices
      {._wkt = "LINESTRING (0 0 0,4 0 4,4 1 4,3 1 3,2 1 2,1 1 1,0 1 0,0 0 0)",
       ._isPlane = true},
      // same outline, one vertex pushed out of the plane
      {._wkt = "LINESTRING (0 0 0,4 0 4,4 1 4,3 1 3,2 1 2.5,1 1 1,0 1 0,0 0 0)",
       ._isPlane = false},
  }};

  for (size_t t = 0; t != test.size(); ++t) {
    std::unique_ptr<Geometry> const g(io::readWkt(test[t]._wkt));
    const LineString               *l = dynamic_cast<LineString *>(g.get());
    BOOST_CHECK_MESSAGE(
        algorithm::isPlane3D<Kernel>(*l, 1.e-9) == test[t]._isPlane,
        std::format("LineString {}: {} {}", t, test[t]._wkt,
                    (test[t]._isPlane ? "is plane" : "isn't plane")));
  }
}

BOOST_AUTO_TEST_CASE(testPlane3DDivideByZeroCrash)
{
  std::unique_ptr<Geometry> degenerate_polygon =
      io::readWkt("POLYGON ((1 -1 -1,1 0.5 0.5,1 0.5 0.5,1 -1 -1))");
  BOOST_CHECK(degenerate_polygon->geometryTypeId() == TYPE_POLYGON);

  // Should return degenerate plane without throwing
  auto degenerate_plane =
      algorithm::plane3D<Kernel>(degenerate_polygon->as<Polygon>());

  // See triangulatePolygon3D for this pattern
  if (algorithm::hasPlane3D<Kernel>(degenerate_polygon->as<Polygon>())) {
    // Should not get here, OR plane3D with Plane3DInexactUnsafe should not
    // divide by zero
    auto div_by_zero_check = algorithm::plane3D<Kernel>(
        degenerate_polygon->as<Polygon>(), algorithm::Plane3DInexactUnsafe());
  }

  std::unique_ptr<Geometry> ok_polygon =
      io::readWkt("POLYGON ((1 0.5 0.5,1.5 1.5 0.5,1.5 0.5 0.5,1 0.5 0.5))");
  BOOST_CHECK(ok_polygon->geometryTypeId() == TYPE_POLYGON);

  BOOST_CHECK(algorithm::hasPlane3D<Kernel>(ok_polygon->as<Polygon>()));

  auto valid_plane = algorithm::plane3D<Kernel>(
      ok_polygon->as<Polygon>(), algorithm::Plane3DInexactUnsafe());
}

BOOST_AUTO_TEST_SUITE_END()

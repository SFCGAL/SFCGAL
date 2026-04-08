// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/intersection.h"
#include "SFCGAL/Geometry.h"
#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/algorithm/equality.h"
#include "SFCGAL/algorithm/intersects.h"
#include "SFCGAL/io/wkt.h"

#include "../../../test_config.h"

#include <boost/test/unit_test.hpp>

using namespace boost::unit_test;
using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(SFCGAL_IntersectionTest)

//
// https://trac.osgeo.org/postgis/ticket/4157
BOOST_AUTO_TEST_CASE(test_postgis_4157)
{
  std::unique_ptr<Geometry> const g1(
      io::readWkt("POLYGON Z (("
                  "122395.299 489126.697 8.61546664325712,"
                  "122389.298 489128.73 8.55588025324629,"
                  "122391.489 489135.198 8.5526708028059,"
                  "122397.49 489133.165 8.61225719281685,"
                  "122395.299 489126.697 8.61546664325712))"));
  std::unique_ptr<Geometry> const g2(
      io::readWkt("POLYHEDRALSURFACE Z ((("
                  "122390.998245685 489133.068537491 0,"
                  "122391.003145022 489133.066423547 0,"
                  "122391.003145022 489133.066423547 10,"
                  "122390.998245685 489133.068537491 10,"
                  "122390.998245685 489133.068537491 0"
                  ")),(("
                  "122391.003145022 489133.066423547 0,"
                  "122383.269575402 489114.842869866 0,"
                  "122383.269575402 489114.842869866 10,"
                  "122391.003145022 489133.066423547 10,"
                  "122391.003145022 489133.066423547 0"
                  ")))"));

  algorithm::intersection3D(*g1, *g2);
}

BOOST_AUTO_TEST_CASE(testIntersectionTINZ)
{
  const std::unique_ptr<Geometry> cylinder = io::readWkt(
      "PolyhedralSurface Z (((150 90 0, 135.35 54.64 0, 100 40 0, 64.64 54.64 "
      "0, 50 90 0, 64.64 125.35 0, 100 140 0, 135.35 125.35 0, 150 90 0)),"
      "((150 90 30, 135.35 125.35 30, 100 140 30, 64.64 125.35 30, 50 90 30, "
      "64.64 54.64 30, 100 40 30, 135.35 54.64 30, 150 90 30)),"
      "((150 90 0, 150 90 30, 135.35 54.64 30, 135.35 54.64 0, 150 90 0)),"
      "((135.35 54.64 0, 135.35 54.64 30, 100 40 30, 100 40 0, 135.35 54.64 "
      "0)),"
      "((100 40 0, 100 40 30, 64.64 54.64 30, 64.64 54.64 0, 100 40 0)),"
      "((64.64 54.64 0, 64.64 54.64 30, 50 90 30, 50 90 0, 64.64 54.64 0)),"
      "((50 90 0, 50 90 30, 64.64 125.35 30, 64.64 125.35 0, 50 90 0)),"
      "((64.64 125.35 0, 64.64 125.35 30, 100 140 30, 100 140 0, 64.64 125.35 "
      "0)),"
      "((100 140 0, 100 140 30, 135.35 125.35 30, 135.35 125.35 0, 100 140 0)),"
      "((135.35 125.35 0, 135.35 125.35 30, 150 90 30, 150 90 0, 135.35 125.35 "
      "0)))");

  const std::unique_ptr<Geometry> cube = io::readWkt(
      "PolyhedralSurface Z (((130 80 0, 80 30 0, 30 80 0, 80 130 0, 130 80 0)),"
      "((130 80 30, 80 130 30, 30 80 30, 80 30 30, 130 80 30)),"
      "((130 80 0, 130 80 30, 80 30 30, 80 30 0, 130 80 0)),"
      "((80 30 0, 80 30 30, 30 80 30, 30 80 0, 80 30 0)),"
      "((30 80 0, 30 80 30, 80 130 30, 80 130 0, 30 80 0)),"
      "((80 130 0, 80 130 30, 130 80 30, 130 80 0, 80 130 0)))");

  BOOST_CHECK(algorithm::intersects3D(*cylinder, *cube));

  std::unique_ptr<Geometry> intersection =
      algorithm::intersection3D(*cylinder, *cube);
  GeometryCollection result = intersection->as<GeometryCollection>();
  BOOST_CHECK(!result.isEmpty());
  BOOST_CHECK_EQUAL(result.numGeometries(), 7);

  const std::unique_ptr<Geometry> expectedGeom1 =
      io::readWkt("LINESTRING Z (57.07 107.07 13.76,57.07 107.07 0.00)");
  BOOST_CHECK(algorithm::almostEqual(result.geometryN(0), *expectedGeom1, 1e2));
}

BOOST_AUTO_TEST_SUITE_END()

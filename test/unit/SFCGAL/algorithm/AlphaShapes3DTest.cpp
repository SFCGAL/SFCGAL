// Copyright (c) 2026-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <boost/format/parsing.hpp>
#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test.hpp>

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/algorithm/alphaShapes3D.h"
#include "SFCGAL/algorithm/isValid.h"
#include "SFCGAL/io/wkt.h"

using namespace SFCGAL;

#include "../../../test_config.h"
// always after CGAL
using namespace boost::unit_test;

BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_AlphaShapes3DTest)

// algorithm::alphaShapes3D

BOOST_AUTO_TEST_CASE(testAlphaShapes3D_Empty)
{
  GeometryCollection emptyCollection;
  emptyCollection.addGeometry(Polygon());
  emptyCollection.addGeometry(Polygon());
  std::unique_ptr<Geometry> emptyAlphaShape3D(
      algorithm::alphaShapes3D(emptyCollection));
  BOOST_CHECK(emptyAlphaShape3D->isEmpty());
}

BOOST_AUTO_TEST_CASE(testAlphaShapes3D_MultiPoint)
{
  // input data
  std::string inputData(SFCGAL_TEST_DIRECTORY);
  inputData += "/data/bunny1000Wkt.txt";
  std::ifstream bunnyFS(inputData.c_str());
  BOOST_REQUIRE(bunnyFS.good());

  std::ostringstream inputWkt;
  inputWkt << bunnyFS.rdbuf();

  std::unique_ptr<Geometry> geomInput(io::readWkt(inputWkt.str()));
  BOOST_REQUIRE(geomInput->is3D());

  // no cavities
  std::string expectedResult(SFCGAL_TEST_DIRECTORY);
  expectedResult += "/data/AlphaShapes3DWkt_expected.txt";
  std::ifstream efsResult(expectedResult.c_str());
  BOOST_REQUIRE(efsResult.good());

  std::string expectedWkt;
  std::getline(efsResult, expectedWkt);

  std::unique_ptr<Geometry> alphaShapes(
      algorithm::alphaShapes3D(geomInput->as<const SFCGAL::Geometry>()));
  BOOST_CHECK_EQUAL(alphaShapes->asText(6), expectedWkt);
}

BOOST_AUTO_TEST_SUITE_END()

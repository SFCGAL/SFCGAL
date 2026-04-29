// Copyright (c) 2026-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test.hpp>

#include "SFCGAL/Geometry.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/algorithm/covers.h"
#include "SFCGAL/algorithm/mergeCoplanarFaces.h"
#include "SFCGAL/io/wkt.h"

using namespace SFCGAL;

#include "../../../test_config.h"
// always after CGAL
using namespace boost::unit_test;

BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_MergeCoplanarFaces)

BOOST_AUTO_TEST_CASE(testMergeCoplanarFaces_Empty)
{
  const PolyhedralSurface emptyGeom;
  BOOST_CHECK(emptyGeom.isEmpty());
  std::unique_ptr<Geometry> result = algorithm::mergeCoplanarFaces(emptyGeom);
  BOOST_CHECK(result->isEmpty());
}

BOOST_AUTO_TEST_CASE(testMergeCoplanarFaces_OneFace)
{
  const std::string faceTINWkt =
      "TIN Z (((0 0 0, 1 0 0, 1 1 0, 0 0 0)), ((0 0 0, 1 1 0, 0 1 0, 0 0 0)))";
  const std::unique_ptr<Geometry> faceTIN = io::readWkt(faceTINWkt);

  std::unique_ptr<Geometry> mergedGeom =
      SFCGAL::algorithm::mergeCoplanarFaces(*faceTIN);
  BOOST_CHECK_EQUAL(mergedGeom->asText(0),
                    "POLYHEDRALSURFACE Z (((1 1 0,0 1 0,0 0 0,1 0 0,1 1 0)))");
}

BOOST_AUTO_TEST_CASE(testMergeCoplanarFaces_Solid)
{
  std::string inputData(SFCGAL_TEST_DIRECTORY);
  inputData += "/data/split3D/houseSolidTooManyFaces.wkt";
  std::ifstream ifs(inputData.c_str());
  BOOST_REQUIRE(ifs.good());

  std::string inputWkt;
  std::getline(ifs, inputWkt);
  const std::unique_ptr<Geometry> houseSolid(io::readWkt(inputWkt));
  BOOST_CHECK(!houseSolid->isEmpty());
  BOOST_CHECK_EQUAL(houseSolid->geometryTypeId(),
                    SFCGAL::GeometryType::TYPE_SOLID);
  BOOST_CHECK_EQUAL(houseSolid->as<Solid>().shellN(0).numPatches(), 9);

  std::unique_ptr<Geometry> houseSimplified =
      SFCGAL::algorithm::mergeCoplanarFaces(*houseSolid);
  BOOST_CHECK_EQUAL(houseSimplified->as<Solid>().shellN(0).numPatches(), 7);

  BOOST_CHECK(algorithm::covers3D(*houseSolid, *houseSimplified));
}

BOOST_AUTO_TEST_SUITE_END()
